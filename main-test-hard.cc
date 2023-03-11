#include "catsc/cats.hh"

#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>
#include <cassert>

/**\file
 * \brief Test on trivial filter
 *
 * This test assumes general correctness of the CATS(C) algorithm: N layers
 * are filled with random "hits". Geometrical filter is not applied, --
 * instead we explicitly imply connections between hits. Algorithm is expected
 * to reveal all the implied chains.
 * */

static void
_usage_info(const char * appName, std::ostream & os) {
    os << "Usage:" << std::endl
       << "   $ "  << appName << " <seed:int> <nLayers:int> <nHits:int> <nMissedLayers:int> <minLength:int>" << std::endl
       ;
}

struct HardHit : public std::unordered_set<const HardHit *> {
    std::pair<size_t, size_t> id;
    HardHit(size_t n1, size_t n2) : id{n1, n2} {}
};
struct HardLayer : public std::unordered_set<const HardHit *> { };
// Custom finder type to derive interface classes
typedef catsc::TrackFinder<const HardHit *> CATSTrackFinder;
// Checks triplet of hits vs "hard" topology
struct HardFilter : public CATSTrackFinder::iTripletFilter {
    bool matches( const HardHit & a
                , const HardHit & b
                , const HardHit & c
                ) const override {
        auto itAB = b.find(&a)
           , itBC = c.find(&b)
           ;
        return itAB != b.end() && itBC != c.end();
    }
};

namespace std {
template <>
struct hash< std::vector<const HardHit *> > {
    std::size_t operator()(const std::vector<const HardHit *> & v) const {
        size_t nComponent = 0;
        size_t r = reinterpret_cast<size_t>(*v.begin());
        for(auto ptr : v) {
            if(!(nComponent++)) continue;
            size_t num = reinterpret_cast<size_t>(ptr);
            num <<= (nComponent%4);
            num >>= (nComponent%4);
            r ^= num;
        }
        return r;
    }
};

std::ostream & operator<<(std::ostream & os, const HardHit & hh) {
    return os << "(" << hh.id.first << "," << hh.id.second
              /*<< ":" << hh.size()*/ << ")";
}

}  // namespace std

// Collector type. Once track found, corresponding flag is raised.
struct Collector : public CATSTrackFinder::iTrackCandidateCollector
                 , public std::unordered_map< std::vector<const HardHit *>, bool >
                 {
protected:
    // recursive helper
    void _fill_from_hit(std::vector<const HardHit *> &, size_t minLength);
public:
    /// gets called vs every found track candidate
    void collect(const cats_HitData_t *, size_t ) override;
    /// used to pre-fill the collector with expected track candidates
    void prepare(const std::vector<HardLayer> &, size_t minLength=3);

    /// Prints state of the collector
    void dump(std::ostream & os) const;
};

void
Collector::_fill_from_hit(std::vector<const HardHit *> & stack, size_t minLength) {
    for(auto toHitPtr : *stack.back()) {
        stack.push_back(toHitPtr);
        if(stack.size() >= minLength && toHitPtr->empty()) emplace(stack, false);
        _fill_from_hit(stack, minLength);
        assert(stack.back() == toHitPtr);
        stack.pop_back();
    }
}

void
Collector::prepare(const std::vector<HardLayer> & layers, size_t minLength) {
    std::vector<const HardHit *> stack;
    for( auto & layer : layers ) {
        for( auto hitPtr : layer ) {
            stack.push_back(hitPtr);
            _fill_from_hit(stack, minLength);
            assert(stack.back() == hitPtr);
            stack.pop_back();
        }
    }
}

void
Collector::dump(std::ostream & os) const {
    for(auto item : *this) {
        os << (item.second ? "1" : "0") << " : ";
        for(auto hitPtr : item.first) {
            os << " " << *hitPtr;
        }
        os << std::endl;
    }
}

void
Collector::collect( const cats_HitData_t * hits_
                  , size_t length
                  ) {
    std::vector<const HardHit *> found;
    for(size_t i = 0; i < length; ++i) {
        const HardHit * hitPtr = reinterpret_cast<const HardHit *>(hits_[i]);
        //std::cout << *hitPtr << std::endl;  // XXX
        found.push_back(hitPtr);
    }
    auto it = find(found);
    bool error = false;
    if(it == end()) {
        error = true;
        std::cerr << "CATS(C) found non-existing candidate: ";
    } else {
        if(it->second) {
            error = true;
            std::cerr << "CATS(C) found same candidate twice: ";
        } else {
            it->second = true;
        }
    }
    if(error) {
        for(auto p : found) {
            std::cerr << *p << ", ";
        }
        std::cerr << std::endl;
    }
}

int
main(int argc, char * argv[]) {
    if(argc != 6) {
        _usage_info(argv[0], std::cerr);
        return 1;
    }
    // set initial parameters
    const int nLayers = std::stoi(argv[2])
            , nHitsPerLayer = std::stoi(argv[3])
            , nMissedLayers = std::stoi(argv[4])
            , minLength = std::stoi(argv[5])
            ;
    // init random number generators
    std::default_random_engine eng{std::stoul(argv[1])};
    std::vector<HardLayer> layers;
    assert(nHitsPerLayer > 1);
    std::uniform_int_distribution<int> urd( nHitsPerLayer*.75
                                          , nHitsPerLayer*1.25
                                          );
    // file to print out evaluation states, as JSON
    FILE * debugJSONStream = fopen("debug.json", "w");
    // Track finder instance to populate
    CATSTrackFinder cats(nLayers, 1000, 10, 10, debugJSONStream);
    // generate "hard" hits topology
    for(int nLayer = 0; nLayer < nLayers; ++nLayer) {
        // number of hits per layer choosen randomly +/- 50% of base value
        int nHitsOnLayer = urd(eng);
        layers.push_back(HardLayer{});
        std::cout << nHitsOnLayer << " hits on layer #" << nLayer << std::endl;
        // populate layer with hits
        for(int nHit = 0; nHit < nHitsOnLayer; ++nHit) {
            HardHit * hit = new HardHit(nLayer, nHit);
            cats.add(nLayer, hit);  // for ptr it does not matter when to insert
            // establish some random backward connections:
            bool isFirst = true;
            for( int toLayer = nLayer - nMissedLayers - 1
               ; toLayer < nLayer
               ; ++toLayer ) {
                if(toLayer < 0) continue;
                std::uniform_int_distribution<int> urd2( 1
                                            , (layers[toLayer].size()-1)/2 );
                int nLinks = urd2(eng);
                std::uniform_int_distribution<int> urd3( 0
                                            , layers[toLayer].size()-1 );
                std::unordered_set<int> appeared;
                for( int nLink = 0; nLink < nLinks; ++nLink ) {
                    auto it = layers[toLayer].begin();
                    int n = urd3(eng);
                    if(appeared.find(n) != appeared.end()) continue;
                    appeared.insert(n);
                    std::advance(it, n);
                    assert(it != layers[toLayer].end());
                    std::cout << " connecting " << *hit << " -> " << **it
                              << " [" << n << "]" << std::endl;
                    if(isFirst) isFirst = false;
                    hit->insert(*it);
                }
            }
            layers.back().insert(hit);
        }
    }
    Collector collector;
    collector.prepare(layers, minLength);

    HardFilter f;
    //fputs("{\"iterations\":[", debugJSONStream);
    cats.collect(f, collector, minLength, nMissedLayers);
    //fputs("]}", debugJSONStream);

    collector.dump(std::cout);  // XXX, initial state

    std::cout << "done." << std::endl;
    return 0;
}
