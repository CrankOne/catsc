#include <fstream>
#include <regex>
#include <iostream>

static const std::regex gRxGraphLine (
    R"~(\s*(\d+)\s*,\s*(\d+)\s*->\s*(\d+)\s*,\s*(\d+)\s*\:\s*(\d+\.\d+)\s*)~");

int
read_graph (std::ifstream & ifs) {
    std::string line;
    size_t nLine = 0;
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
        std::cout << "(" << fromL << "," << fromN << ") -> ("
                  << toL << "," << toN << ") : " << w << std::endl;
    }
    return 0;
}

int
main(int argc, char * argv[]) {
    std::ifstream ifs(argv[1]);
    read_graph(ifs);
    return -0;
}
