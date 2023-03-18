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
#include <cstring>

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
struct HardFilterUnweighted : public CATSTrackFinder::iTripletFilter {
    bool matches( const HardHit & a
                , const HardHit & b
                , const HardHit & c
                ) const override {
        auto itAB = b.find(&a)
           , itBC = c.find(&b)
           ;
        //std::cout << "  testing a=" << a << ", b=" << b << ", c=" << c    // XXX
        //          << " => " << (itAB != b.end()) << (itBC != c.end())       // XXX
        //          << std::endl;                                             // XXX
        return itAB != b.end() && itBC != c.end();
    }
};
// Checks triplet of hits vs "hard" topology
struct HardFilterWeighted : public CATSTrackFinder::iWeightedTripletFilter {
    cats_Weight_t weight( const HardHit & a
                        , const HardHit & b
                        , const HardHit & c
                        ) const override {
        auto itAB = b.find(&a)
           , itBC = c.find(&b)
           ;
        if( itAB == b.end() || itBC == c.end() ) return 0;
        //std::cout << "  weighting a=" << a << ", b=" << b << ", c=" << c    // XXX
        //          << " => " << itAB->second * itBC->second                  // XXX
        //          << std::endl;                                             // XXX
        return itAB->second * itBC->second;
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

// A (nasty) type for test cases, to be filled from file. Represents a sequence
// of found candidates of the hierarchical form:
//  #found/(is-found:bool, hits:[(layer#, hit#)])
typedef
    std::vector<std::pair<bool, std::vector<std::pair<size_t, size_t>>>>
    ExpectedConnections;

// Collector type. Once track found, corresponding flag is raised.
struct Collector : public CATSTrackFinder::iTrackCandidateCollector
                 , public std::unordered_map< std::vector<const HardHit *>, bool >
                 {
public:
    bool verbose;
    ExpectedConnections expected;
    bool failed;
protected:
    size_t _nCollected;
    // recursive helper
    void _fill_from_hit(std::vector<const HardHit *> &, size_t minLength);
public:
    Collector() : verbose(false), failed(false), _nCollected(0) {}
    ~Collector();
    /// gets called vs every found track candidate
    void collect(const cats_HitData_t *, size_t ) override;
    /// used to pre-fill the collector with expected track candidates
    void prepare(const std::vector<HardLayer> &, size_t minLength=3);

    /// Prints state of the collector
    void dump(std::ostream & os) const;

    void reset() {
        _nCollected = 0;
        failed = false;
        for(auto & entry : *this) entry.second = false;
        expected.clear();
    }
};

Collector::~Collector() {
}

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
    if(verbose) std::cout << "  track candidate:";
    for(size_t i = 0; i < length; ++i) {
        const HardHit * hitPtr = reinterpret_cast<const HardHit *>(hits_[i]);
        if(verbose) std::cout << " " << *hitPtr;
        found.push_back(hitPtr);
    }
    if(verbose) std::cout << std::endl;

    bool error = false;
    if(!empty()) {
        auto it = find(found);
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

    if(!expected.empty()) {
        if(_nCollected >= expected.size()) {
            std::cerr << "Comparison with expected ignored (candidate number"
                " exceeds expected)." << std::endl;
            ++_nCollected;
            failed = true;
            return;
        }
        bool errorOnEntry = false;
        auto & entry = expected[_nCollected];
        if(entry.second.size() != found.size()) {
            std::cerr << "size mismatch on #" << _nCollected << ":" << std::endl;
            errorOnEntry = true;
        }
        for(size_t i = 0; i < entry.second.size(); ++i) {
            if( found[i]->id != entry.second[i] ) {
                std::cerr << i << "-th entry mismatch on #" << _nCollected << ":" << std::endl;
                errorOnEntry = true;
                break;
            }
        }
        if(errorOnEntry) {
            std::cerr << "  expected: ";
            for(const auto & hit : entry.second ) {
                std::cerr << hit.first << "," << hit.second << " ";
            }
            std::cerr << std::endl;
            std::cerr << "     found: ";
            for(const HardHit * hit : found ) {
                std::cerr << hit->id.first << "," << hit->id.second << " ";
            }
            std::cerr << std::endl;
            failed = true;
        } else {
            assert(!entry.first);  // found same candidate twice
            entry.first = true;
        }
    }
    ++_nCollected;
}

//                          * * *   * * *   * * *

static std::vector<HardLayer>
_fill_random( std::vector<HardLayer> & layers //CATSTrackFinder & cats
            , unsigned long seed
            , int nLayers
            , int nHitsPerLayer
            , int nMissedLayers
            , int minLength
            , bool verbose
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
        if(verbose)
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
                    //std::cout << " connecting " << *hit << " -> " << **it  // XXX
                    //          << " [" << w << "]" << std::endl;            // XXX
                    if(isFirst) isFirst = false;
                    hit->emplace(*it, w);
                }
            }
            //layers.back().push_back(hit);
        }
    }
    return layers;
}

//                          * * *   * * *   * * *

typedef
    std::tuple<std::string, cats_LayerNo_t, cats_LayerNo_t, bool>  // strategy, length, missed, weighted
    TestCaseKey;

namespace std {
template <> struct hash<TestCaseKey> {
    std::size_t operator()(const TestCaseKey & k) const {
      return ((std::hash<std::string>()(std::get<0>(k))
               ^ (std::get<1>(k) << 1)) << 1)
               ^ (std::get<2>(k) << 8)
               ^ (std::get<3>(k) ? 0x24 : 0x188)
               ;
    }
};
}  // namespace std
class ExpectedCases : public std::unordered_map< TestCaseKey
                                               , ExpectedConnections
                                               > {};

static const std::regex gRxGraphLine (
    R"~(\s*(\d+)\s*,\s*(\d+)\s*->\s*(\d+)\s*,\s*(\d+)\s*\:\s*(\d+\.\d+)\s*)~");
static const std::regex gRxStartCaseLine (
    R"~(\s*CASE\s+([a-z]+)\s+length\s+(\d+)\s+missed\s+(\d+)\s+(un)?weighted\s*)~");
static const std::regex gRxCaseItem (
    R"~(\s*(\d+)\s*,\s*(\d+)\s*)~");

static size_t
_read_test_cases ( std::ifstream & ifs
                 , std::vector<HardLayer> & layers_
                 , ExpectedCases & expConns
                 ) {
    std::string line;
    size_t nLine = 0;
    std::map<size_t, std::map<size_t, HardHit *>> layers;
    // this is for case reading
    std::string cStrategy;
    cats_LayerNo_t nMissed = 0, minLength = 0;
    bool weighted = false;

    size_t nTestCaseItem;
    ExpectedConnections expConnsStrat;  // reentrant fill buffer
    // iterate over lines
    while (std::getline (ifs, line)) {
        ++nLine;
        line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), line.end());
        if (line.empty()) continue;  // omit empty lines
        if ('#' == line[0]) continue;  // omit comment lines
        std::smatch match;
        if(!cStrategy.empty()) {  // if we're in test case block, read line as item
            if( line == "ENDCASE" ) {  // note: sensitive to whitespaces
                auto ir = expConns.emplace( TestCaseKey{
                                              cStrategy
                                            , nMissed
                                            , minLength
                                            , weighted
                                            }
                                    , expConnsStrat );
                if(!ir.second) {
                    char errbf[128];
                    snprintf(errbf, sizeof(errbf), "Test case on line %zu"
                                " repeats one of the above.", nLine);
                    throw std::runtime_error(errbf);
                }
                cStrategy.clear();
                nMissed = minLength = 0;
                weighted = false;
                continue;
            }
            std::smatch sm;
            std::string::const_iterator sit(line.begin());
            std::vector<std::pair<size_t, size_t>> items;
            while(std::regex_search(sit, line.cend(), sm, gRxCaseItem)) {
                std::vector<std::string> stoks(sm.begin(), sm.end());
                std::pair<size_t, size_t> item(std::stoi(stoks[1]), std::stoi(stoks[2]));
                items.push_back(item);
                sit = sm.suffix().first;
            }
            assert(nTestCaseItem == expConnsStrat.size());
            expConnsStrat.push_back(decltype(expConnsStrat)::value_type{false, items});
            ++ nTestCaseItem;
            continue;
        }
        if(std::regex_match(line, match, gRxStartCaseLine)) {  // start test case block
            std::vector<std::string> toks (match.begin(), match.end());
            cStrategy = toks [1];
            minLength = std::stoi(toks[2]);
            nMissed = std::stoi(toks[3]);
            weighted = toks[4].empty();
            nTestCaseItem = 0;
            expConnsStrat.clear();
            continue;
        }
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
        //std::cout << " connecting " << *fromHitPtr << " -> " << *toHitPtr
        //          << " [" << w << "]" << std::endl;
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
    bool weighted;
    FILE * debugJSONStream;
    bool verbose;
    bool noBruteforce;

    struct {
        // this is used for randomly-generated mode
        struct {
            long seed;
            int nLayers;
            int nHitsPerLayerBase;
        } rndGenOpts;
        struct {
            std::string filename;
            bool testAll;
        } readOpts;
    } mtdDep;
};

static void
_set_defaults( AppCfg & cfg ) {
    cfg.mode = '\0';
    cfg.nMissedLayers = 1;
    cfg.nMinLength = 3;
    cfg.weighted = false;
    cfg.debugJSONStream = NULL;
    cfg.verbose = false;
    cfg.noBruteforce = false;

    cfg.mtdDep.rndGenOpts.seed = 1337;
    cfg.mtdDep.rndGenOpts.nLayers = 5;
    cfg.mtdDep.rndGenOpts.nHitsPerLayerBase = 7;

    cfg.mtdDep.readOpts.filename = "";
    cfg.mtdDep.readOpts.testAll = false;
}

static void
_usage_info( std::ostream & os
           , const char * appName
           , AppCfg & cfg
           ) {
    os << "CATS(C) test application." << std::endl
       << "Usage:" << std::endl
       << "    (1) $ " << appName << " -i <filename> -A [-v]" << std::endl
       << "    (2) $ " << appName << " -i <filename> [COMMON_OPTS]" << std::endl
       << "    (3) $ " << appName << " -s <seed> -N <#layers> -n <#hits-per-layer> [COMMON_OPTS]" << std::endl
       << "In modes (1) and (2) will read topology from file provided as argument to"
          " `-i'. In mode (3) will generate random topology respecting seed"
          " given by `-s' with average number of hits `-n' on `-N' layers."
          " In all cases app will perform then finding procedure for the topology"
          " (read or generated). Th strategy and"
          " constrains are expected to be provided as [COMMON_OPTS]."
          " When running in mode (2) with [COMMON_OPTS] matching one of the"
          " test cases, the collected candidates will be automatically tested"
          " against found test case (if any). In mode"
          " (1) expected output will be tested against all"
          " test cases provided with the given file (so COMMON_OPTS are"
          " not needed)." << std::endl
       << "COMMON_OPTS are:" << std::endl
       << "  -m <method>, required option, defines collection method. Must"
          " be one of: \"exessive\", \"strict\", \"weighted-strict\","
          " \"longest\", \"weighted-longest\"."
       << std::endl
       << "  -l <min-length=" << cfg.nMinLength << "> restriction on minimal"
          " number of hits on the found track candidates." << std::endl
       << "  -w <#hits=" << cfg.nMissedLayers << "> number of tolerated missed"
          " (inefficient) layers." << std::endl
       << "  -W enables weighted connection procedure." << std::endl
       << "  -D <debugJSONStream> when set, debug JSON will be generated by"
          " given path." << std::endl
       << "  -c turns off brute-force collection procedure for read or built"
          " graph. Useful for profiling (as brute-forcing graph dominates the"
          " the CATS."
       << std::endl
       ;
}

static int
_configure_app( int argc, char * argv[]
              , AppCfg & cfg
              ) {
    int opt;
    while ((opt = getopt(argc, argv, "hvacWi:s:N:n:m:l:w:D:")) != -1) {
        switch (opt) {
        case 'h':
            _usage_info(std::cout, argv[0], cfg);
            return 1;
        // common opts
        case 'v':
            cfg.verbose = true;
            break;
        case 'm':
            cfg.method = optarg;
            break;
        case 'l':
            cfg.nMinLength = std::stoi(optarg);
            break;
        case 'w':
            cfg.nMissedLayers = std::stoi(optarg);
            break;
        case 'W':
            cfg.weighted = true;
            break;
        case 'c':
            cfg.noBruteforce = true;
            break;
        case 'D':
            cfg.debugJSONStream = fopen(optarg, "w");
            break;
        #if 0
        case 't':
            //cfg.mtdDep.readOpts.testCases.push_back(optarg);
            if(!strcmp(optarg, "all")) {
                cfg.mtdDep.readOpts.testCases.push_back(TestCaseKey{"all", 0, 0});
            } else {
                // parse test case description
                std::regex tcRx(R"~(([a-z]+)-(\d+)-(\d+))~");
                std::smatch m;
                std::string expr(optarg);
                if(!std::regex_match(expr, m, tcRx)) {
                    std::cerr << "Test reference from cmd-line \"" << optarg << "\""
                        << " doesn't match the pattern." << std::endl;
                    return 1;
                }
                std::vector<std::string> toks(m.begin(), m.end());
                cfg.mtdDep.readOpts.testCases.push_back(
                       TestCaseKey(toks[1], std::stoi(toks[2]), std::stoi(toks[3])) );
            }
            break;
        #endif
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
        case 'a':
            if('\0' == cfg.mode) cfg.mode = 'f';
            if(cfg.mode != 'f') {
                std::cerr << "App config error: conflicting keys"
                    " (contradicting modes)." << std::endl;
                return -2;
            }
            cfg.mtdDep.readOpts.testAll = true;
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
    ExpectedCases expConns;

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
        if(cfg.verbose)
            std::cout << "Reading \"" << cfg.mtdDep.readOpts.filename << "\":" << std::endl;
        size_t totalHits = _read_test_cases(ifs, layers, expConns);
        if(cfg.verbose) {
            std::cout << "Read " << layers.size() << " layers with "
                      << totalHits << " hits total." << std::endl;
            std::cout << "Test cases discovered:" << std::endl;
            for(auto & testItem : expConns) {
                std::cout << " * \"" << std::get<0>(testItem.first) << "\", missed="
                    << std::get<1>(testItem.first) << ", length="
                    << std::get<2>(testItem.first) << std::endl;
            }
        }
    } else if(cfg.mode == 'r') {  // generate random topology
        _fill_random( layers
                    , cfg.mtdDep.rndGenOpts.seed
                    , cfg.mtdDep.rndGenOpts.nLayers
                    , cfg.mtdDep.rndGenOpts.nHitsPerLayerBase
                    , cfg.nMissedLayers
                    , cfg.nMissedLayers
                    , cfg.verbose
                    );
    } else {
        std::cerr << "App error: no mode (nothing to do)." << std::endl;
        return 1;
    }

    if(cfg.verbose) {
        std::cout << "Hits on layers:" << std::endl;
        size_t n = 0;
        for( auto layerEntry : layers ) {
            std::cout << "   #" << n++ << " of " << layerEntry.size() << " hits"
                << std::endl;
        }
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
    collector.verbose = cfg.verbose;
    if(!cfg.noBruteforce)
        collector.prepare(layers);

    std::unordered_set<TestCaseKey> keys;
    if(!cfg.mtdDep.readOpts.testAll) {
        keys.emplace(TestCaseKey{cfg.method, cfg.nMissedLayers, cfg.nMinLength, cfg.weighted});
    } else {
        if(expConns.empty()) {
            std::cerr << "WARNING: can't perform all tests as file does not"
                " contain any test cases." << std::endl;
        }
        std::transform( expConns.begin(), expConns.end()
                      , std::inserter(keys, keys.end())
                      , [](auto & p) {return p.first;}
                      );
    }
    bool hadError = false
       , testPerformed = false
       ;
    for(auto k : keys) {
        collector.reset();
        auto control = expConns.find(k);
        if(control == expConns.end()) {
            if(cfg.verbose)
                std::cout << "No matching test case (" << std::get<0>(k)
                    << ", length=" << std::get<2>(k) << ", missed="
                    << std::get<1>(k) << ")."
                    << std::endl;
        } else {
            if(cfg.verbose)
                std::cout << "Considering test case (" << std::get<0>(k)
                    << ", length=" << std::get<2>(k) << ", missed="
                    << std::get<1>(k) << ") of "
                    << control->second.size() << " entries."
                    << std::endl;
            collector.expected = control->second;
        }
        bool evaluated = false;
        const cats_LayerNo_t thisNMissing = std::get<1>(k)
                           , thisLength = std::get<2>(k)
                           ;
        const std::string thisMethod = std::get<0>(k);
        const bool weighted = std::get<3>(k);

        //HardFilterUnweighted filter;
        CATSTrackFinder::iTripletFilterBase * filterBase = nullptr;
        if(!weighted) {
            filterBase = new HardFilterUnweighted();
        } else {
            filterBase = new HardFilterWeighted();
        }

        try {
            if(!weighted)
                evaluated = cats.evaluate(*static_cast<HardFilterUnweighted *>(filterBase), thisNMissing);
            else
                evaluated = cats.evaluate(*static_cast<HardFilterWeighted *>(filterBase), thisNMissing);
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

        if( thisMethod == "excessive" ) {
            cats.collect_excessive(collector, thisLength);
        } else if( thisMethod == "moderate" ) {
            cats.collect(collector, thisLength);
        } else if( thisMethod == "strict" ) {
            cats.collect_strict(collector, thisLength);
        } else if(thisMethod == "longest") {
            cats.collect_longest(collector, thisLength);
        } else if(thisMethod == "winning") {
            cats.collect_winning(collector, thisLength);
        } else {
            std::cerr << "Unknown collection"
                " method: \"" << thisMethod << "\"." << std::endl;
            return 1;
        }

        if(cfg.verbose)
            collector.dump(std::cout);

        hadError |= collector.failed;
        if(!collector.expected.empty()) {
            testPerformed = true;
            size_t nEntry = 0;
            for(const auto & entry : collector.expected) {
                if(entry.first) {
                    ++nEntry;
                    continue;
                }
                hadError = true;
                std::cerr << "expected entry #" << nEntry << " missed:";
                for(auto p : entry.second) {
                    std::cerr << " " << p.first << "," << p.second;
                }
                std::cerr << std::endl;
                ++nEntry;
            }
        }
        if(filterBase)
            delete filterBase;
    }  // modes loop
    if(cfg.debugJSONStream) {
        fclose(cfg.debugJSONStream);
        cfg.debugJSONStream = nullptr;
    }

    if(cfg.verbose) {
        if(hadError)
            std::cerr << "There were errors." << std::endl;
        else if(testPerformed)
            std::cout << "Test ok." << std::endl;
        else
            std::cout << "No test case foreseen." << std::endl;
    }

    for( auto layerEntry : layers ) {
        for( auto item : layerEntry ) {
            delete item;
        }
    }

    return hadError ? 1 : 0;
}

