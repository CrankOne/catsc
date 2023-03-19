# Generic C/C++ implementation of cellular automata for track finding

This code is based on the algorithm of cellular automaton evolution
proposed for track finding by I.Abt et al [1]. Algorithm is designed to
perform initial track finding for experiments in High energy physics, yet
its topological properties may imply a broader usage.

This implementation does does not assume any particular geometry and uses
abstract (yet highly efficient) leveled graph to exploit various strategies of
[Depth-first search](https://en.wikipedia.org/wiki/Depth-first_search) to yield
what is assumed to be track candidates.

# Usage and library design

All core functions are written in pure C to gain competetive performance for
various low-level applications (like high-level trigger). Certain
optimizations are applied to reduce heap usage and gain performance. User code
pursuing performance should rely on pure C interface.

For less demanding applications a bit more laconic and expressive C++ interface
based on dynamic polymorphism is provided. Yet, expect some performance penalty
due to virtual functions resolution.

For all cases user code has to define following routines:

- A triplet *filter* -- a routine providing answer on whether three *hits*
  instances can be a part of the *track candidate*.
- A *collector* -- a routine or function that shall be invoked on every found
  track candidate.

A *filter* can return discrete answer or a weight (as a positive floating point
number). In first case insertion order deterministically affects the collecting
order.

There also different strategies to iterate over track candidates (in order of
decreased redundancy):

- **excessive** strategy will visit all possible track candidates, including
  their sub-sequences till the minimal length is reached (longest come first).
  It won't visit same candidate twice and is more efficient than naive
  brute-force combinatorics.
- **moderate** is similar to **excessive** but won't visit same triplets unless
  they're concurrent (at same depth)
- **strict** strategy will provide candidates with intersections, but omit exact
  sub-sequences
- **longest** resolves concurrency in favour of longest branch. Especially
  useful in combination with weighted filter.
- **winning** resolves concurrency exclusively. This strategy effectively,
  combined with weighted filter is equivalent to one proposed in the original
  article. Provides candidates without repeating *pairs* of hits.

Note, that due to potential redundancy CATS itself does not keep or return
found track candidates. It is reasonably cheap, however to re-iterate them
from user code.

## C++ interface

```cpp
// It is convenient to typedef the main CATS class for further usage; this
// template must be parameterised with user hit type
typedef catsc::TrackFinder<UserHitType> TrackFinder;

// Inherit and implement interface for filter
class Filter : public TrackFinder::iTripletFilter {
    bool matches (const UserHitType & a, const UserHitType & b, const UserHitType & c) const override {
        // ... return whether a,b,c is a valid triplet
    }
};
// Alternatively, implement weighted filter
//class UserFilter : public TrackFinder::iWeightedTripletFilter {
//    cats_Weight_t weight (const UserHitType & a, const UserHitType & b, const UserHitType & c) const override {
//        // ... return weight of a,b,c triplet; return a value <=0 if triplet
//        //     can not be considered as valid at all
//    }
//};

// Implement track candidate collector
class Collector : public TrackFinder::iTrackCandidateCollector {
    void collect (const cats_HitData_t * hits, size_t nHits) override {
        // for sake of efficiency a pointer to void (void*) is exposed here,
        // so one needs a cast to get actual hits data:
        for(size_t i = 0; i < nHits; ++i) {
            const UserHit * hitPtr = reinterpret_cast<const HardHit *>(hits[i]);
            // ... copy hits or perform additional validation stages to reduce
            // output
        }
    }
};

// Use above as follows:

// create main reentrant track finding instance parameterised with
// layers number:
TrackFinder cats(nLayers);
// Instantiate user fitter and collector reentrant instances
Filter filter;
Collector collector;

while (eventsAvailable) {  // within the phys.events loop
    while (hitsAvailableInEvent) {
        // populate layers with hits data
        cats.add (layerNo, hit);
    }
    // evaluate collector instance, parameterised with filter (weighted or not)
    // and number of tolerable missing layers
    cats.evaluate(filter, nMissedLayers);
    // collect the track candidates with collectro insance using favoured
    // strategy of required min length
    cats.collect_strict(collector, minLength);
    // reset automaton before next event
    cats.reset();
}
```

## C interface

```c
// implement unweighted filter function
int user_filter_func (cats_HitData_t a_, cats_HitData_t b_, cats_HitData_t c_, void * userdata) {
    // shall return 1 if triplet is valid and 0 otherwice. Cast to user hit
    // data type ptr to access actual data
    struct UserHitType * a = (const struct UserHitType *) a_;
    // ...
}
// or use weighted function of signature (returns floating point)
//cats_Weight_t user_filter_func_w (cats_HitData_t a_, cats_HitData_t b_, cats_HitData_t c_, void * userdata) {
//...
//}

// implement collecting function
void user_collecting_func (const cats_HitData_t * tc, size_t tcLen, void * userdata) {
    for(size_t tcN = 0; tcN < tcLen; ++tcN) {
        struct UserHitType * a = *((const struct UserHitType * const *) tc)[tcN];
        // ...
    }
}

// initialize automaton with maximum layers number with
cats_Layers * cats = cats_layers_create(6);
// create workspace for automaton (allocation pools)
cats_CellsPool * pool = cats_cells_pool_create(4);

while (eventsAvailable) {  // within the phys.events loop
    while (hitsAvailableInEvent) {
        // populate layers with pointers to hit data
        cats_layer_add_point (cats, nLayer, &hit);
    }
    // initialize graph with your filter function; provide pointer user data,
    // if need. `nMissedLayers` is max number of missing layer
    // (inefficient/missed detectors)
    cats_connect (cats, pool, user_filter_func, &user_filter_data, nMissedLayers);
    // ...or use weighted connection function if weighted filter is in use
    //cats_connect_w (cats, pool, user_filter_func_w, &user_filter_data, nMissedLayers);
    // evaluate automaton. 2-nd parameter is the debug stream, set to NULL if
    // you don't need one
    cats_evaluate (cats, NULL);
    // collect the tracks with the visiter object using certain visiting
    // strategy.
    cats_visit_dfs_excessive (cats, minLength, user_collecting_func, &userCollectingData );
    // ...or use weighted visitor if weighted filter function was used:
    //cats_visit_dfs_excessive_w (cats, minLength, user_collecting_func, &userCollectingData );
    // reset graph and pool before re-using same workspace
    cats_cells_pool_reset (_layers, _cells, _softLimitCells);
    cats_layers_reset (_layers, _softLimitHits);
}
// delete workspaces
cats_cells_pool_delete(pool);
cats_layers_delete(cats);
```

Refer to `tests/main-gpl.cc` for example of usage of plain C API (yet,
incomplete in terms of reentrant usage).

## Build

Standard CMake build procedure for ou-of-source build:

    $ cd build
    $ cmake ..
    $ make -j4 install

Note, that user might be interested in customizing type for hit identification
and/or floating point type to describe spatial coordinates -- for that purpose
once can set CMake variables `CATS_HIT_ID_TYPE` and `CATS_COORDINATE_TYPE`
with `-D...=...` respictively.

## CMake project integration

For CMake-based projects:

    find_package(catsc)
    # ...
    target_link_libraries(${na64std_LIB} PUBLIC catsc::catsc )

For other types of build system a `pkg-config` file is exported by given user
prefix.

# TODO list

* Doublet filter
* User-driven hits lookup
* C++ API based on static polymorphism
* Parallel CPU implementation
* OpenCV implementation
* Rough math library (approximate trigonometry and exponentiation)

# References

[1] Abt, I.; Emeliyanov, D.; Kisel, I.; Masciocchi, S.; CATS: a cellular
automaton for tracking in silicon for the HERA-B vertex detector // Nuclear
Instruments and Methods in Physics Research A 489 (2002) 389–405

