#include "catsc/cats.hh"

#include <stdexcept>
#include <algorithm>

//
// Main C++ wrapper class

CATSTrackFinder::CATSTrackFinder( cats_LayerNo_t nLayers
                                , size_t softLimitCells
                                , size_t softLimitHits
                                , size_t softLimitRefs
                                )
    : _nLayers(nLayers)
    , _softLimitCells(softLimitCells)
    , _softLimitHits(softLimitHits)
    , _softLimitRefs(softLimitRefs)
    {
    if( nLayers < 3 ) {
        throw std::runtime_error("Bad number of layers requested (<3).");
    }
    _layers = cats_layers_create(nLayers);
    if(!_layers) throw std::bad_alloc();
    _cells = cats_cells_pool_create(nLayers);
    if(!_cells) throw std::bad_alloc();
}

CATSTrackFinder::~CATSTrackFinder() {
    cats_layers_delete(_layers);
    cats_cells_pool_delete(_cells);
}

void
CATSTrackFinder::add( cats_LayerNo_t nLayer, cats_HitID_t id
                    , cats_Float_t x, cats_Float_t y, cats_Float_t z) {
    int rc = cats_layer_add_point(_layers, nLayer, x, y, z, id);
    if(rc) throw std::bad_alloc();
}

static int
_filter_f_wrapper( cats_HitID_t id1, const cats_Float_t * c1
                 , cats_HitID_t id2, const cats_Float_t * c2
                 , cats_HitID_t id3, const cats_Float_t * c3
                 , void * filter_ ) {
    // NOTE: it appears that resolving virtual ptr on vtable here drains
    //       performance to significant extent.
    return reinterpret_cast<CATSTrackFinder::iTripletFilter*>(filter_)
        ->matches( id1, c1, id2, c2, id3, c3) ? 1 : 0;
}

static void
_collect_f_wrapper( cats_HitID_t * hitIDs, size_t nHits, void * collector_ ) {
    reinterpret_cast<CATSTrackFinder::iTrackCandidateCollector*>(collector_)
        ->collect(hitIDs, nHits);
}

void
CATSTrackFinder::collect( iTripletFilter & filter
                        , iTrackCandidateCollector & collector
                        , unsigned int minLength
                        , unsigned int nMissingLayers
                        ) {
    if( cats_evolve( _layers, _cells
                   , _filter_f_wrapper
                   , &filter
                   , nMissingLayers ) ) {
        throw std::bad_alloc();
    }
    cats_for_each_track_candidate(_layers, minLength, _collect_f_wrapper, &collector);
    cats_cells_pool_reset(_layers, _cells, _softLimitCells);
}

void
CATSTrackFinder::reset() {
    cats_layers_reset(_layers, _softLimitHits);
}

//
// An utility class to provide "longest unique" sequences

void
LongestUniqueTrackCollector::collect(cats_HitID_t * hitIDs, size_t nHitIDs) {
    std::set<cats_HitID_t> newSet(hitIDs, hitIDs + nHitIDs);
    std::set<cats_HitID_t> diff;
    bool hasNewHits = true;
    for( const auto & cSet : _collected ) {
        // if given set is a subset of any of the already considered ones,
        // omit it
        std::set_difference( newSet.begin(), newSet.end()
                           , cSet.begin(), cSet.end()
                           , std::inserter(diff, diff.begin())
                           );
        if( ! diff.empty() ) {
            // new element didn't bring new hits with respect to one of
            // the already collected
            return;
        }
    }
    // if we are here, the hits sequence is unique and shall be considered
    // as a track candidate
    _collected.push_back(newSet);
    _consider_track_candidate(hitIDs, nHitIDs);
}

