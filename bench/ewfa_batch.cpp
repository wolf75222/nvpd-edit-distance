// Baseline path-producing batch runners for the eWFA-GPU comparison.
#include "chps/cuda.hpp"
#include "edlib.h"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

int main(int argc, char** argv) try {
    if (argc != 3) throw std::invalid_argument("ewfa_batch {edlib|nvpd-cuda} pairs.seq");
    const std::string algorithm = argv[1];
    if (algorithm != "edlib" && algorithm != "nvpd-cuda") throw std::invalid_argument("unknown algorithm");
    std::ifstream in(argv[2]);
    std::vector<std::pair<std::string,std::string>> pairs;
    std::string a, b;
    while (std::getline(in,a)) {
        if (!std::getline(in,b) || !a.starts_with('>') || !b.starts_with('<'))
            throw std::invalid_argument("invalid pair");
        pairs.emplace_back(a.substr(1),b.substr(1));
    }
    if (pairs.empty()) throw std::invalid_argument("empty or unreadable input");
    const auto start = std::chrono::steady_clock::now();
    std::vector<int> scores;
    {
        chps::CudaWorkspace workspace;
        for (const auto& [query,target] : pairs) {
            if (algorithm == "edlib") {
                auto result = edlibAlign(query.data(),query.size(),target.data(),target.size(),
                                         edlibNewAlignConfig(-1,EDLIB_MODE_NW,EDLIB_TASK_PATH,nullptr,0));
                if (result.status != EDLIB_STATUS_OK) throw std::runtime_error("edlib failed");
                scores.push_back(result.editDistance);
                edlibFreeAlignResult(result);
            } else {
                scores.push_back(workspace.distance(query,target,true).distance);
            }
        }
    }
    const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout << "seconds " << seconds << '\n';
    for (int d : scores) std::cout << "distance " << d << '\n';
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}
