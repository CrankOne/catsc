#include "catsc/cats.hh"

#include <stdexcept>
#include <algorithm>

namespace catsc {

//
// Main C++ wrapper class

void
BaseTrackFinder::c_f_wrapper_collect( cats_HitID_t * hitIDs, size_t nHits, void * collector_ ) {
    reinterpret_cast<BaseTrackFinder::iTrackCandidateCollector*>(collector_)
        ->collect(hitIDs, nHits);
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

}  // namespace catsc

