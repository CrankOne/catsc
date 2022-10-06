#pragma once

#include "catsc/cats.h"

#include <vector>
#include <set>

/**\brief A C++ CATS track finding algorithm representation.
 *
 * This class wraps the C routines into C++ class providing enclosed lifecycle
 * of reentrant objects:
 *
 * 1. Track finder instance has to be initialized with certain number of layers
 * 2. Hits (space points) have to be added with `add()` method
 * 3. Track candidates with certain filter cand collector has to be collected,
 *    possibly not once (it is possible to apply different filters)
 * 4. Hits can be re-set and 2,3 can be repeated.
 * 5. Destroyed track finder instance frees all the memory in reentrant structs.
 *
 * Cells topology is defined by filter.
 *
 * The collector will be supplied with all the track candidates found for
 * certain topology. Note, that track candidates can be intersecting meaning
 * that for sequences (track candidates) of more than three hits (or whatever
 * is specified to `collect()`'s `minLength`) a sub-sequences will appear. To
 * de-duplicate, consider using of `LongestUniqueTrackCollector` helper class.
 * */
class CATSTrackFinder {
public:
    /// Base class for hit triplet filtering functor
    struct iTripletFilter {
        /// Should return whether triplet is possibly describes a segment of
        /// a track based on IDs and coordinates of three points. Lifetime of
        /// the provided ptrs is valid until `reset()` call.
        virtual bool matches( cats_HitID_t, const cats_Float_t *
                            , cats_HitID_t, const cats_Float_t *
                            , cats_HitID_t, const cats_Float_t *
                            ) const = 0;
    };

    /// Base class for track candidate collecting functors
    struct iTrackCandidateCollector {
        /// Shall consider/collect the track candidate defined as ordered
        /// array of IDs of given size. Lifetime of provided array is not
        /// persistant.
        virtual void collect(cats_HitID_t *, size_t) = 0;
    };
private:
    /// Number of layers this instance was initialized to handle
    cats_LayerNo_t _nLayers;
    /// Ptr to C-struct representing the layers structure
    cats_Layers * _layers;
    /// Ptr to pool of C-structs representing connection between hits ("cells")
    cats_CellsPool * _cells;
    size_t _softLimitCells, _softLimitHits, _softLimitRefs;
public:
    /// Constructs track finder instance for certain number of "layers"
    CATSTrackFinder( cats_LayerNo_t nLayers
                   , size_t softLimitCells=10000
                   , size_t softLimitHits=100
                   , size_t softLimitRefs=1000
                   );
    /// dtr, pretty strightforward one.
    ~CATSTrackFinder();

    /// Adds point (a hit) to be considered on certain layer
    void add( cats_LayerNo_t nLayer, cats_HitID_t id
            , cats_Float_t, cats_Float_t, cats_Float_t);

    /// Returns number of layers the finder was initialized to handle
    cats_LayerNo_t n_layers() const { return _nLayers; }

    /// Builds the cells relying on provided filter, evaluates the algorithm
    /// on cells, and iterates over the all set of found tracks using
    /// collector.
    ///
    /// \note Added hits are erased once all the track candidates are collected.
    void collect( iTripletFilter &
                , iTrackCandidateCollector &
                , unsigned int minLength
                , unsigned int nMissingLayers
                );
    /// Drops associated hits
    void reset();
};

/**\brief A collector shim filtering only longest unique hit sequences
 *
 * This shimming class shall provide only the sequences of hits that are not
 * part of each other. It exploits the property of CATS that at first the
 * collector is provided with longest sequences and then their sub-seequences
 * go. If sub-sequence is fully contained by one of the previously appeared
 * track candidate(s), it will be omitted.
 *
 * In comparison with original DFS search proposed by CATS papers it shall
 * correctly identify the cases with doubling hits.
 * */
class LongestUniqueTrackCollector : public CATSTrackFinder::iTrackCandidateCollector {
private:
    ///\brief A storage of considered unique track candidates
    std::vector<std::set<cats_HitID_t>> _collected;
protected:
    ///\brief This method will be called only with unique track candidates
    virtual void _consider_track_candidate(cats_HitID_t * hitIDs, size_t nHitIDs) = 0;
public:
    ///\brief Implements lookup for duplicates
    virtual void collect(cats_HitID_t * hitIDs, size_t nHitIDs) override;
    /// Shall be called at the end of each event to drop the cache with track
    /// candidates being already considered
    virtual void reset() {
        _collected.clear();
    }
};

