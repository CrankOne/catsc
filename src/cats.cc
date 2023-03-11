#include "catsc/cats.hh"

#include <stdexcept>
#include <algorithm>

#include <iostream>  // XXX

namespace catsc {

//
// Main C++ wrapper class

void
BaseTrackFinder::c_f_wrapper_collect( const cats_HitData_t * hitIDs, size_t nHits, void * collector_ ) {
    reinterpret_cast<BaseTrackFinder::iTrackCandidateCollector*>(collector_)
        ->collect(hitIDs, nHits);
}

//
// An utility class to provide "longest unique" sequences

void
LongestUniqueTrackCollector::collect(const cats_HitData_t * hits, size_t nHitIDs) {
    std::set<cats_HitData_t> newSet(hits, hits + nHitIDs);
    std::set<cats_HitData_t> diff;
    for( const auto & cSet : _collected ) {
        // if given set is a subset of any of the already considered ones,
        // omit it
        std::set_difference( newSet.begin(), newSet.end()
                           , cSet.begin(), cSet.end()
                           , std::inserter(diff, diff.begin())
                           );
        if( diff.empty() ) {
            // new element didn't bring new hits with respect to one of
            // the already collected
            std::cout << "  ... track declined" << std::endl;  // XXX
            return;
        }
        diff.clear();
    }
    // if we are here, the hits sequence is unique and shall be considered
    // as a track candidate
    _collected.push_back(newSet);
    _consider_track_candidate(hits, nHitIDs);
}

}  // namespace catsc

