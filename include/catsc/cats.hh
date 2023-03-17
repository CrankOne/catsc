#pragma once

#include "catsc/cats.h"

#include <vector>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <cassert>

#include <iostream>  // XXX

namespace catsc {

/// Not for direct use; provides C function to be used in template and
/// collector interface base
struct BaseTrackFinder {
    /// Base class for track candidate collecting functors
    struct iTrackCandidateCollector {
        /// Shall consider/collect the track candidate defined as ordered
        /// array of IDs of given size. Lifetime of provided array is not
        /// persistant.
        virtual void collect(const cats_HitData_t *, size_t) = 0;
        /// Called after all the tracks has been collected.
        virtual void done() {}
    };

    /// C callback implem to forward hits into collector
    static void c_f_wrapper_collect( const cats_HitData_t *, size_t nHits, void * );
};

template<typename DataT> struct HitDataTraits {
    typedef std::vector<DataT *> CacheType;
    typedef const DataT & HitInfo_t;
    static DataT * get_data_ptr(CacheType & cacheRef, const DataT & obj) {
        DataT * localCopyPtr = new DataT(obj);
        cacheRef.push_back(localCopyPtr);
        return localCopyPtr;
    }
    static void reset(CacheType & cacheRef) { cacheRef.clear(); }
    typedef const DataT * HitInfoPtr_t;
};

template<typename DataT> struct HitDataTraits<DataT*> {
    struct CacheType {};
    typedef const DataT & HitInfo_t;
    static DataT * get_data_ptr(CacheType &, DataT * ptr) { return ptr; }
    static void reset(CacheType & cacheRef) {}
    typedef const DataT * HitInfoPtr_t;
};

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
 *
 * "Soft limit" is the number of certain entities allowed in the event and may
 * affects the performance in some way. The idea is that reentrant structures
 * usually do not free their memory at `reset()` call to decrease the number
 * of actual reallocations. However, practically it is somitemes possible that
 * event will cause major increase of such pools. To soften an impact of this
 * kind of allocations, a limit is set for these entities. If pool exceeds this
 * limit it will be actually `free()`'d, otherwise the memory is
 * kept (reserved).
 * */
template<typename HitDataT>
class TrackFinder : protected HitDataTraits<HitDataT>::CacheType
                  , public BaseTrackFinder
                  {
public:
    struct iTripletFilterBase {
        const bool weighted;
        iTripletFilterBase(bool isWeighted) : weighted(isWeighted) {}
        virtual ~iTripletFilterBase() {}
    };
    /// Base class for hit triplet filtering functor
    struct iTripletFilter : public iTripletFilterBase {
        iTripletFilter() : iTripletFilterBase(false) {}
        /// Should return whether triplet is possibly describes a segment of
        /// a track based on IDs and coordinates of three points. Lifetime of
        /// the provided ptrs is valid until `reset()` call.
        virtual bool matches( typename HitDataTraits<HitDataT>::HitInfo_t
                            , typename HitDataTraits<HitDataT>::HitInfo_t
                            , typename HitDataTraits<HitDataT>::HitInfo_t
                            ) const = 0;
    };
    /// Base class for weighted hit triplet functor
    struct iWeightedTripletFilter : public iTripletFilterBase {
        iWeightedTripletFilter() : iTripletFilterBase(true) {}
        /// Sould return weight value for given triplet
        virtual float weight( typename HitDataTraits<HitDataT>::HitInfo_t
                            , typename HitDataTraits<HitDataT>::HitInfo_t
                            , typename HitDataTraits<HitDataT>::HitInfo_t
                            ) const = 0;
    };
private:
    /// Number of layers this instance was initialized to handle
    cats_LayerNo_t _nLayers;
    /// Ptr to C-struct representing the layers structure
    cats_Layers * _layers;
    /// Ptr to pool of C-structs representing connection between hits ("cells")
    cats_CellsPool * _cells;
    /// Soft limits for allocations
    size_t _softLimitCells, _softLimitHits, _softLimitRefs;

    /// C callback implem applying the filter
    static int c_f_wrapper_filter( cats_HitData_t c1
                                 , cats_HitData_t c2
                                 , cats_HitData_t c3
                                 , void * filter_
                                 );
    /// C callback implem applying weighted filter
    static cats_Weight_t c_f_wrapper_wfilter( cats_HitData_t c1
                                            , cats_HitData_t c2
                                            , cats_HitData_t c3
                                            , void * filter_
                                            );

    /// If set, inctance is ready to produce tracks
    bool _evaluated;
    /// Set to number of missing layers for which the automaton was evaluated
    cats_LayerNo_t _lastNMissing;
    /// If set to non-null pointer, states for every iteration will be printed
    /// as JSON into the file
    FILE * _debugJSONStream;
    /// Internal flag, used to print debug info when stream ptr is given
    bool _isFirstHit;
    /// Internal flag set after `collect()` call; used to re-set internal flags
    /// before possible repeatative `collect()` calls
    bool _wasCollected;
    /// Last return code of C procedure
    int _rc;
protected:
    /// Connects the layers and evaluates the automaton.
    ///
    /// \throw `std::bad_alloc` on memory error
    bool _evaluate(iTripletFilterBase & filter_, cats_LayerNo_t nMissingLayers) {
        if( _debugJSONStream ) fputs("},\"its\":[", _debugJSONStream);

        if(!filter_.weighted) {
            auto & filter = static_cast<iTripletFilter&>(filter_);
            _rc = cats_connect( _layers
                              , _cells
                              , TrackFinder<HitDataT>::c_f_wrapper_filter
                              , &filter
                              , nMissingLayers
                              );
        } else {
            auto & filter = static_cast<iWeightedTripletFilter&>(filter_);
            _rc = cats_connect_w( _layers
                                , _cells
                                , TrackFinder<HitDataT>::c_f_wrapper_wfilter
                                , &filter
                                , nMissingLayers
                                );
        }

        if(_rc) {
            if(_rc == CATSC_RC_EMPTY_GRAPH ) {
                _evaluated = true;
                return false;  // filter discriminated all
            }
            if(!(_rc & CATSC_ERROR_RUNTIME_LOGIC)) {
                throw std::bad_alloc();
            } else {
                char errBf[64];
                snprintf( errBf, sizeof(errBf), "CATS(C) runtime error: %#x."
                        , _rc);
                throw std::runtime_error(errBf);
                // ^^^ TODO: verbose descriptions/dedicated exceptions for
                //     error codes
            }
        }
        if( 0 != cats_evaluate( _layers, _debugJSONStream )) {
            if(_debugJSONStream) fputs("]}", _debugJSONStream);
            _evaluated = false;
            return false;
        }
        if(_debugJSONStream) fputs("]}", _debugJSONStream);
        _evaluated = true;
        return true;
    }

    ///\brief Must be called before repeatative `collect()` calls
    void _reset_collection_flags()
        { reset_collection_flags(_cells); }
public:
    /// Constructs track finder instance for certain number of "layers"
    ///
    /// \throw `std::bad_alloc` if failed to allocate memory for layers or cells
    TrackFinder( cats_LayerNo_t nLayers
               , size_t softLimitCells=10000
               , size_t softLimitHits=100
               , size_t softLimitRefs=1000
               , FILE * debugJSONStream=nullptr
               ) : _nLayers(nLayers)
                 , _softLimitCells(softLimitCells)
                 , _softLimitHits(softLimitHits)
                 , _softLimitRefs(softLimitRefs)
                 , _evaluated(false)
                 , _lastNMissing(0)
                 , _debugJSONStream(debugJSONStream)
                 , _isFirstHit(true)
                 , _rc(0)
                 {
        if( nLayers < 3 ) {
            throw std::runtime_error("Bad number of layers requested (<3).");
        }
        _layers = cats_layers_create(nLayers);
        if(!_layers) throw std::bad_alloc();  // failed to allocate layers
        _cells = cats_cells_pool_create(nLayers);
        if(!_cells) throw std::bad_alloc();  // failed to allocate cells pool
        if( _debugJSONStream ) fputs("{\"dict\":{", _debugJSONStream);
    }
    /// dtr, pretty strightforward one.
    ~TrackFinder() {
        cats_layers_delete(_layers);
        cats_cells_pool_delete(_cells);
    }

    const int c_retcode() const { return _rc; }

    /// Adds point (a hit) to be considered on certain layer
    ///
    /// \throw `std::bad_alloc` if could not allocate memory for hit
    void add( cats_LayerNo_t nLayer, HitDataT data ) {
        if( _evaluated )
            throw std::runtime_error("Adding hits to instance with"
                    " inapropriate state (after CATS evaluated).");  // missing reset() call?
        auto p = HitDataTraits<HitDataT>::get_data_ptr(*this, data);
        assert(p);
        int rc = cats_layer_add_point(_layers, nLayer, p);
        if(rc != 0) throw std::bad_alloc();
        if(_debugJSONStream) {
            if(!_isFirstHit) fputc(',', _debugJSONStream); else _isFirstHit = false;
            int nHit = cats_layer_n_points(_layers, nLayer);
            fprintf( _debugJSONStream
                   , "\"%p\":\"(%d,%d)\""
                   , p, nLayer, nHit - 1 );
        }
    }

    /// Returns number of layers the finder was initialized to handle
    cats_LayerNo_t n_layers() const { return _nLayers; }

    /// Returns number of hits added to layer N
    size_t n_points(cats_LayerNo_t nLayer) { return cats_layer_n_points(_layers, nLayer); }

    /// Evaluates automaton; prepares connection graph for tracks collection
    bool evaluate( iTripletFilterBase & filter
                 , cats_LayerNo_t nMissingLayers ) {
        if(_evaluated) {
            cats_cells_pool_reset( _layers
                                 , _cells
                                 , _softLimitCells
                                 );
            _evaluated = false;
        }
        _lastNMissing = nMissingLayers;
        return _evaluate(filter, _lastNMissing);
    }

    ///\brief Collects all track candidates permitted by the filter, including
    ///       all subsequences.
    ///
    /// Forwards execution to `cats_for_each_track_candidate_excessive()` that
    /// may generate highly excessive result (see docs).
    void collect_excessive( iTrackCandidateCollector & collector
                          , cats_LayerNo_t minLength
                          ) {
        if(!_evaluated)
            throw std::runtime_error("CATS was not evaluated, unable to collect.");
        if(_wasCollected)
            _reset_collection_flags();
        int rc = cats_visit_dfs_excessive_w(_layers, minLength  // TODO: not weighted!
                , c_f_wrapper_collect, &collector);
        collector.done();
        if(rc) throw std::runtime_error("collect() error");  // TODO: elaborate
        _wasCollected = true;
    }

    ///\brief Collects all track candidates permitted by the filter
    ///
    /// Forwards execution to `cats_for_each_track_candidate()` that might
    /// be excessive (see docs).
    void collect( iTrackCandidateCollector & collector
                , cats_LayerNo_t minLength
                ) {
        if(!_evaluated)
            throw std::runtime_error("CATS was not evaluated, unable to collect.");
        if(_wasCollected)
            _reset_collection_flags();
        int rc = cats_visit_dfs_moderate(_layers, minLength/*, _lastNMissing*/
                , c_f_wrapper_collect, &collector);
        if(rc) throw std::runtime_error("collect() error");  // TODO: elaborate
        collector.done();
    }

    ///\brief Collects all track candidates permitted by the filter
    ///
    /// Forwards execution to `cats_for_each_track_candidate_strict()` that might
    /// be insufficient or excessive (see docs).
    void collect_strict( iTrackCandidateCollector & collector
                       , cats_LayerNo_t minLength
                       ) {
        if(!_evaluated)
            throw std::runtime_error("CATS was not evaluated, unable to collect.");
        if(_wasCollected)
            _reset_collection_flags();
        int rc = cats_visit_dfs_strict(_layers, minLength/*, _lastNMissing*/
                , c_f_wrapper_collect, &collector);
        if(rc) throw std::runtime_error("collect() error");  // TODO: elaborate
        collector.done();
    }

    ///\brief Collects longest track candidates permitted by the filter
    ///
    /// Forwards execution to `cats_for_each_longest_track_candidate()` that
    /// might be insufficient or imprecise (see docs).
    void collect_longest( iTrackCandidateCollector & collector
                        , cats_LayerNo_t minLength
                        ) {
        if(!_evaluated)
            throw std::runtime_error("CATS was not evaluated, unable to collect.");
        if(_wasCollected)
            _reset_collection_flags();
        int rc = cats_visit_dfs_longest( _layers, minLength, _lastNMissing
                                       , c_f_wrapper_collect, &collector);
        if(rc) throw std::runtime_error("collect() error");  // TODO: elaborate
        collector.done();
    }

    ///\brief Collects winning track candidates permitted by the filter
    ///
    /// Forwards execution to `cats_for_each_winning_track_candidate()` that
    /// might be insufficient or imprecise (see docs).
    void collect_winning( iTrackCandidateCollector & collector
                        , cats_LayerNo_t minLength
                        ) {
        if(!_evaluated)
            throw std::runtime_error("CATS was not evaluated, unable to collect.");
        if(_wasCollected)
            _reset_collection_flags();
        cats_visit_dfs_winning( _layers, minLength, _lastNMissing
                              , c_f_wrapper_collect, &collector);
        collector.done();
    }

    /// Drops associated hits
    void reset() {
        cats_cells_pool_reset(_layers, _cells, _softLimitCells);
        _isFirstHit = true;
        cats_layers_reset(_layers, _softLimitHits);
        HitDataTraits<HitDataT>::reset(*this);
        _evaluated = false;
        cats_cells_pool_reset(_layers, _cells, _softLimitCells);
    }
};

template<typename HitDataT> int
TrackFinder<HitDataT>::c_f_wrapper_filter( cats_HitData_t c1
                                         , cats_HitData_t c2
                                         , cats_HitData_t c3
                                         , void * filter_
                                         ) {
    return reinterpret_cast<typename TrackFinder<HitDataT>::iTripletFilter *>(filter_)
        ->matches( *reinterpret_cast<typename HitDataTraits<HitDataT>::HitInfoPtr_t>(c1)
                 , *reinterpret_cast<typename HitDataTraits<HitDataT>::HitInfoPtr_t>(c2)
                 , *reinterpret_cast<typename HitDataTraits<HitDataT>::HitInfoPtr_t>(c3)
                 ) ? 1 : 0;
}

template<typename HitDataT> cats_Weight_t
TrackFinder<HitDataT>::c_f_wrapper_wfilter( cats_HitData_t c1
                                          , cats_HitData_t c2
                                          , cats_HitData_t c3
                                          , void * filter_
                                          ) {
    return reinterpret_cast<typename TrackFinder<HitDataT>::iWeightedTripletFilter *>(filter_)
        ->weight( *reinterpret_cast<typename HitDataTraits<HitDataT>::HitInfoPtr_t>(c1)
                , *reinterpret_cast<typename HitDataTraits<HitDataT>::HitInfoPtr_t>(c2)
                , *reinterpret_cast<typename HitDataTraits<HitDataT>::HitInfoPtr_t>(c3)
                );
}

#if 0
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
 *
 * \todo Remove, deprecated in favor of different strategies
 * */
class LongestUniqueTrackCollector : public BaseTrackFinder::iTrackCandidateCollector {
private:
    ///\brief A storage of considered unique track candidates, used for fast lookup
    std::vector<std::set<cats_HitData_t>> _collected;
    ///\brief Storage of currently longest unique track candidates
    std::vector<std::vector<cats_HitData_t>> _collectedData;
protected:
    ///\brief This method will be called only with unique track candidates
    virtual void _consider_track_candidate(const cats_HitData_t * hitIDs, size_t nHitIDs) = 0;
public:
    ///\brief Implements lookup for duplicates
    void collect(const cats_HitData_t * hitIDs, size_t nHitIDs) override;
    /// Runs `_consider_track_candidate()` on selected longest candidates
    void done() override;
    /// Shall be called at the end of each event to drop the cache with track
    /// candidates being already considered
    virtual void reset() {
        _collected.clear();
        _collectedData.clear();
    }
};
#endif

}  // namespace catsc

