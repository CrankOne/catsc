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
    assert(!newSet.empty());
    for( auto it = _collected.begin(); it != _collected.end(); ++it ) {
        const auto & cSet = *it;
        assert(!cSet.empty());
        if(cSet.size() > newSet.size()) {
            // size of new set is less or equal to current --
            // try to locate every element of new set in current. In
            // case new set is a subset of existing, all the "new" elements
            // must be found at existing and we shall ignore new
            bool isSubset = true;
            for(auto newEl : newSet) {
                if(cSet.find(newEl) != cSet.end()) continue;
                // new set brings at least one unique element -- skip further
                // check
                isSubset = false;
                break;
            }
            if(isSubset) return;  // new is a subset, ignore
        } else if( cSet.size() < newSet.size() ) {
            // current set is smaller than new -- check if current is a subset
            // and substitute current for new
            bool isSubset = true;
            for(auto curEl : cSet) {
                if( newSet.find(curEl) != newSet.end() ) continue;
                isSubset = false;
                break;
            }
            if(isSubset) {
                // current set is subset of new; impose new set and that's it
                std::swap(*it, newSet);
                return;
            }
        }
        #ifndef NDEBUG
        else {
            // Algorithm guarantees that we can not have identical sets of the
            // same size visited twice.
            std::set<cats_HitData_t> diff;
            std::set_symmetric_difference( newSet.begin(), newSet.end()
                                         , cSet.begin(), cSet.end()
                                         , std::inserter(diff, diff.begin())
                                         );
            assert(!diff.empty());
        }
        #endif
    }
    // if we are here, the hits sequence is unique and shall be considered
    // as a track candidate
    _collected.push_back(newSet);
    _consider_track_candidate(hits, nHitIDs);
}

}  // namespace catsc

