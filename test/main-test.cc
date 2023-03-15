#include "catsc/cats.hh"

#include <getopt.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <map>

// Geometry-agnostic hit repr -- stores only hit label and weighted connections
struct HardHit : public std::unordered_map<const HardHit *, float> {
    std::pair<size_t, size_t> id;
    HardHit(size_t n1, size_t n2) : id{n1, n2} {}
};
namespace std {
std::ostream & operator<<(std::ostream & os, const HardHit & hh) {
    return os << "(" << hh.id.first << "," << hh.id.second
              /*<< ":" << hh.size()*/ << ")";
}
}  // namespace std
// Geometry-agnostic layer repr -- stores only hits
struct HardLayer : public std::vector<const HardHit *> { };
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
        std::cout << "  testing a=" << a << ", b=" << b << ", c=" << c    // XXX
                  << " => " << (itAB != b.end()) << (itBC != c.end())       // XXX
                  << std::endl;                                             // XXX
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
        stack.push_back(toHitPtr.first);
        if(stack.size() >= minLength && toHitPtr.first->empty()) emplace(stack, false);
        _fill_from_hit(stack, minLength);
        assert(stack.back() == toHitPtr.first);
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
    std::cout << "  track candidate:";
    for(size_t i = 0; i < length; ++i) {
        const HardHit * hitPtr = reinterpret_cast<const HardHit *>(hits_[i]);
        std::cout << " " << *hitPtr;
        found.push_back(hitPtr);
    }
    std::cout << std::endl;
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

//                          * * *   * * *   * * *

static std::vector<HardLayer>
_fill_random( std::vector<HardLayer> & layers //CATSTrackFinder & cats
            , unsigned long seed
            , int nLayers
            , int nHitsPerLayer
            , int nMissedLayers
            , int minLength
            ) {
    // init random number generators
    std::default_random_engine eng{seed};
    assert(nHitsPerLayer > 1);
    std::uniform_int_distribution<int> urd( nHitsPerLayer*.5
                                          , nHitsPerLayer*1.5
                                          );
    // generate "hard" hits topology
    for(int nLayer = 0; nLayer < nLayers; ++nLayer) {
        // number of hits per layer choosen randomly +/- 50% of base value
        int nHitsOnLayer = urd(eng);
        layers.push_back(HardLayer{});
        std::cout << nHitsOnLayer << " hits on layer #" << nLayer << std::endl;
        // populate layer with hits
        for(int nHit = 0; nHit < nHitsOnLayer; ++nHit) {
            HardHit * hit = new HardHit(nLayer, nHit);
            layers.back().push_back(hit);  // for ptr it does not matter when to insert
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
                std::uniform_real_distribution<double> urd4f( 0, 1 );
                std::unordered_set<int> appeared;
                for( int nLink = 0; nLink < nLinks; ++nLink ) {
                    auto it = layers[toLayer].begin();
                    int n = urd3(eng);
                    if(appeared.find(n) != appeared.end()) continue;
                    appeared.insert(n);
                    std::advance(it, n);
                    assert(it != layers[toLayer].end());
                    double w = urd4f(eng);
                    std::cout << " connecting " << *hit << " -> " << **it
                              << " [" << w << "]" << std::endl;
                    if(isFirst) isFirst = false;
                    hit->emplace(*it, w);
                }
            }
            layers.back().push_back(hit);
        }
    }
    return layers;
}

//                          * * *   * * *   * * *

static const std::regex gRxGraphLine (
    R"~(\s*(\d+)\s*,\s*(\d+)\s*->\s*(\d+)\s*,\s*(\d+)\s*\:\s*(\d+\.\d+)\s*)~");

static size_t
_read_graph (std::ifstream & ifs, std::vector<HardLayer> & layers_) {
    std::string line;
    size_t nLine = 0;
    std::map<size_t, std::map<size_t, HardHit *>> layers;
    while (std::getline (ifs, line)) {
        ++nLine;
        if (line.empty()) continue;
        if ('#' == line[0]) continue;
        std::smatch match;
        if(!std::regex_match(line, match, gRxGraphLine)) {
            std::cerr << "line #" << nLine << " invalid" << std::endl;
            continue;
        }
        std::vector<std::string> toks (match.begin(), match.end());
        int fromL = std::stoi (toks[1])
          , fromN = std::stoi (toks[2])
          , toL = std::stoi (toks[3])
          , toN = std::stoi (toks[4])
          ;
        float w = std::stof (toks[5]);

        //std::cout << "    (" << fromL << "," << fromN << ") -> ("
        //          << toL << "," << toN << ") : " << w << std::endl;

        HardHit * fromHitPtr = nullptr; {
            auto fromLayerIt = layers.emplace(fromL, decltype(layers)::mapped_type());
            auto fromHitIt = fromLayerIt.first->second.find(fromN);
            if(fromHitIt == fromLayerIt.first->second.end()) {
                fromHitPtr = new HardHit(fromL, fromN);
                fromLayerIt.first->second.emplace(fromN, fromHitPtr);
            } else {
                fromHitPtr = fromHitIt->second;
            }
        }
        
        HardHit * toHitPtr = nullptr;
        {
            auto toLayerIt = layers.emplace(toL, decltype(layers)::mapped_type());
            auto toHitIt = toLayerIt.first->second.find(toN);
            if(toHitIt == toLayerIt.first->second.end()) {
                toHitPtr = new HardHit(toL, toN);
                toLayerIt.first->second.emplace(toN, toHitPtr);
            } else {
                toHitPtr = toHitIt->second;
            }
        }
        assert(fromHitPtr);
        assert(toHitPtr);
        //fromHitPtr->emplace(toHitPtr, w);
        toHitPtr->emplace(fromHitPtr, w);  // ???
        std::cout << " connecting " << *fromHitPtr << " -> " << *toHitPtr
                  << " [" << w << "]" << std::endl;
    }
    layers_.reserve(layers.rbegin()->first + 1);
    size_t orderlyLayer = 0;
    size_t nTotalHits = 0;
    for(auto layerEntry : layers) {
        if(orderlyLayer != layers_.size()) {
            char errBf[128];
            snprintf( errBf, sizeof(errBf)
                    , "Error filling layer #%zu: probably"
                    " layer(s) of previous orderly number(s) was/were missed."
                    , orderlyLayer );
            throw std::runtime_error(errBf);
        }
        HardLayer layerCopy;
        layerCopy.reserve(layerEntry.second.rbegin()->first + 1);
        size_t orderlyNHit = 0;
        for(auto hitEntry : layerEntry.second) {
            if(hitEntry.first != layerCopy.size()) {
                char errBf[128];
                snprintf( errBf, sizeof(errBf)
                        , "Error filling layer #%zu with hit %zu: probably"
                        " hit(s) of previous orderly number(s) was/were missed."
                        , orderlyLayer, orderlyNHit );
                throw std::runtime_error(errBf);
            }
            layerCopy.push_back(hitEntry.second);
            ++orderlyNHit;
            ++nTotalHits;
        }
        layers_.push_back(layerCopy);
        ++orderlyLayer;
    }

    return nTotalHits;
}

//                          * * *   * * *   * * *

struct AppCfg {
    char mode;
    std::string method;
    int nMissedLayers
      , nMinLength
      ;
    FILE * debugJSONStream;

    struct {
        // this is used for randomly-generated mode
        struct {
            long seed;
            int nLayers;
            int nHitsPerLayerBase;
        } rndGenOpts;
        struct {
            std::string filename;
        } readOpts;
    } mtdDep;
};

static void
_set_defaults( AppCfg & cfg ) {
    cfg.mode = '\0';
    cfg.nMissedLayers = 1;
    cfg.nMinLength = 3;
    cfg.debugJSONStream = NULL;

    cfg.mtdDep.rndGenOpts.seed = 1337;
    cfg.mtdDep.rndGenOpts.nLayers = 5;
    cfg.mtdDep.rndGenOpts.nHitsPerLayerBase = 7;

    cfg.mtdDep.readOpts.filename = "";
}

static void
_usage_info( std::ostream & os
           , const char * appName
           , AppCfg & cfg
           ) {
    os << "CATS(C) test application." << std::endl
       << "Usage:" << std::endl
       << "    (1) $ " << appName << " -i <filename> [COMMON_OPTS]" << std::endl
       << "    (2) $ " << appName << " -s <seed> -N <#layers> -n <#hits-per-layer> [COMMON_OPTS]" << std::endl
       << "In mode (1) will read topology from file provided as argument to"
          " `-i'. In mode (2) will generate random topology respecting seed"
          " given by `-s' with average number of hits `-n' on `-N' layers."
          " In both cases app will perform finding procedure for the topology"
          " (read or generated)." << std::endl
       << "COMMON_OPTS are:" << std::endl
       << "  -m <method>, required option, defines collection method. Must"
          " be one of: \"exessive\", \"strict\", \"weighted-strict\","
          " \"longest\", \"weighted-longest\"."
       << std::endl
       << "  -l <min-length=" << cfg.nMinLength << "> restriction on minimal"
          " number of hits on the found track candidates." << std::endl
       << "  -w <#hits=" << cfg.nMissedLayers << "> number of tolerated missed"
          " (inefficient) layers." << std::endl
       << "  -D <debugJSONStream> when set, debug JSON will be generated by"
          " given path."
       << std::endl
       ;
}

static int
_configure_app( int argc, char * argv[]
              , AppCfg & cfg
              ) {
    int opt;
    while ((opt = getopt(argc, argv, "hi:s:N:n:m:l:w:D:")) != -1) {
        switch (opt) {
        case 'h':
            _usage_info(std::cout, argv[0], cfg);
            return 1;
        // common opts
        case 'm':
            cfg.method = optarg;
            break;
        case 'l':
            cfg.nMinLength = std::stoi(optarg);
            break;
        case 'w':
            cfg.nMissedLayers = std::stoi(optarg);
            break;
        case 'D':
            cfg.debugJSONStream = fopen(optarg, "w");
            break;
        // input file opts
        case 'i':
            if('\0' == cfg.mode) cfg.mode = 'f';
            if(cfg.mode != 'f') {
                std::cerr << "App config error: conflicting keys"
                    " (contradicting modes)." << std::endl;
                return -2;
            }
            cfg.mtdDep.readOpts.filename = optarg;
            break;
        // random gen
        case 's':
            if('\0' == cfg.mode) cfg.mode = 'r';
            if(cfg.mode != 'r') {
                std::cerr << "App config error: conflicting keys"
                    " (contradicting modes)." << std::endl;
                return -3;
            }
            cfg.mtdDep.rndGenOpts.seed = std::stol(optarg);
            break;
        case 'N':
            if('\0' == cfg.mode) cfg.mode = 'r';
            if(cfg.mode != 'r') {
                std::cerr << "App config error: conflicting keys"
                    " (contradicting modes)." << std::endl;
                return -3;
            }
            cfg.mtdDep.rndGenOpts.nLayers = std::stoi(optarg);
            break;
        case 'n':
            if('\0' == cfg.mode) cfg.mode = 'r';
            if(cfg.mode != 'r') {
                std::cerr << "App config error: conflicting keys"
                    " (contradicting modes)." << std::endl;
                return -3;
            }
            cfg.mtdDep.rndGenOpts.nHitsPerLayerBase = std::stoi(optarg);
            break;
        default: /* '?' */
            std::cerr << "App config error: unknown option." << std::endl;
            return -1;
        }
    }
    if(cfg.mode == '\0') {
        std::cerr << "App config error: can't figure out mode (what to do)"
            " by given command line args." << std::endl;
        return -4;
    }
    return 0;
}

int
main(int argc, char * argv[]) {
    AppCfg cfg;
    _set_defaults(cfg);
    int rc = _configure_app(argc, argv, cfg);
    if(rc < 0) return 1;
    if(rc > 0) return 0;

    std::vector<HardLayer> layers;

    if(cfg.mode == 'f') {  // read topology from file
        if(cfg.mtdDep.readOpts.filename.empty()) {
            std::cerr << "File name not set." << std::endl;
            return 1;
        }
        std::ifstream ifs(cfg.mtdDep.readOpts.filename);
        if(!ifs) {
            std::cerr << "Bad file path: \""
                      << cfg.mtdDep.readOpts.filename << "\""
                      << std::endl;
            return 1;
        }
        std::cout << "Reading \"" << cfg.mtdDep.readOpts.filename << "\":" << std::endl;
        size_t totalHits = _read_graph(ifs, layers);
        std::cout << "Read " << layers.size() << " layers with "
                  << totalHits << " hits total." << std::endl;
    } else if(cfg.mode == 'r') {  // generate random topology
        _fill_random( layers
                    , cfg.mtdDep.rndGenOpts.seed
                    , cfg.mtdDep.rndGenOpts.nLayers
                    , cfg.mtdDep.rndGenOpts.nHitsPerLayerBase
                    , cfg.nMissedLayers
                    , cfg.nMissedLayers );
    } else {
        std::cerr << "App error: no mode (nothing to do)." << std::endl;
        return 1;
    }

    std::cout << "Hits on layers:" << std::endl;
    size_t n = 0;
    for( auto layerEntry : layers ) {
        std::cout << "   #" << n++ << " of " << layerEntry.size() << " hits"
            << std::endl;
    }

    CATSTrackFinder cats(layers.size(), 10000, 100, 1000, cfg.debugJSONStream);
    // fill automaton
    size_t nl = 0;
    for( auto layerEntry : layers ) {
        for( auto item : layerEntry ) {
            cats.add( nl, item );
        }
        ++nl;
    }
    Collector collector;
    collector.prepare(layers);

    HardFilter filter;
    bool evaluated = false;
    try {
        evaluated = cats.evaluate(filter, cfg.nMissedLayers);
    } catch( std::bad_alloc & e ) {
        std::cerr << "Bad allocation: " << e.what()
                  << "; last C return code: " << cats.c_retcode()
                  << std::endl;
        return 1;
    }
    if(!evaluated) {
        std::cerr << "CATS was not evaluated (empty graph)." << std::endl;
        return 1;
    }

    if( cfg.method == "excessive" ) {
        cats.collect_excessive(collector, cfg.nMinLength);
    } else if( cfg.method == "moderate" ) {
        cats.collect(collector, cfg.nMinLength);
    } else if( cfg.method == "strict" ) {
        cats.collect_strict(collector, cfg.nMinLength);
    } else if(cfg.method == "longest") {
        cats.collect_longest(collector, cfg.nMinLength);
    } else if(cfg.method == "winning") {
        cats.collect_winning(collector, cfg.nMinLength);
    } else {
        std::cerr << "Unknown collection"
            " method: \"" << cfg.method << "\"." << std::endl;
        return 1;
    }

    //std::cout << "Collected " << ... << " of " << collector.size() << " candidates:"
    //    << std::endl;
    collector.dump(std::cout);

    std::cout << "Test done." << std::endl;

    if(cfg.debugJSONStream) {
        fclose(cfg.debugJSONStream);
        cfg.debugJSONStream = nullptr;
    }

    return 0;
}

