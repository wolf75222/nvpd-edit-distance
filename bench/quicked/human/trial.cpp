#include "astarpa.h"
#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/quicked.h"
#include <chrono>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

struct Alignment {
    int score;
    std::string cigar, source;
    int retries = 0;
};

Alignment align(std::string_view method, const std::string& a, const std::string& b) {
    if (method == "biwfa") {
        auto result = chps::biwfa(a,b,true);
        return {result.distance,std::move(result.cigar),"biwfa"};
    }
    if (method == "quicked-certified") {
        auto result = chps_quicked_certified(a.data(), a.size(), b.data(), b.size());
        return {result.distance,{},result.source,result.retries};
    }
    if (method != "astarpa2-full") throw std::invalid_argument("unknown method");
    uint8_t* raw = nullptr;
    uintptr_t length = 0;
    const auto score = astarpa2_full(reinterpret_cast<const uint8_t*>(a.data()), a.size(),
                                     reinterpret_cast<const uint8_t*>(b.data()), b.size(), &raw, &length);
    std::unique_ptr<uint8_t, decltype(&astarpa_free_cigar)> owner(raw,astarpa_free_cigar);
    if (!raw || score > INT_MAX) throw std::runtime_error("invalid A*PA2 result");
    return {static_cast<int>(score),std::string(reinterpret_cast<const char*>(raw),length),"astarpa2-full"};
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const auto pair = chps::read_fasta(argv[1]);
    if (pair.size() != 2) return 2;
    const std::string& a = pair[0].bases;
    const std::string& b = pair[1].bases;
    if (std::max(a.size(), b.size()) > INT_MAX/100 || a.empty() || b.empty() ||
        a.find_first_not_of("ACGT")!=a.npos || b.find_first_not_of("ACGT")!=b.npos) return 2;
    const int warm = align(argv[2],a,b).score;
    const auto start = std::chrono::steady_clock::now();
    const auto result = align(argv[2],a,b);
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now()-start;
    if (warm < 0 || result.score != warm) return 1;
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status,line)) {
        if (line.rfind("VmHWM:",0)==0) {
            const auto rss = std::stol(line.substr(6));
            if (rss <= 0) return 1;
            std::cout << std::setprecision(12) << elapsed.count() << ' ' << rss << ' ' << result.score << ' '
                      << result.source << ' ' << result.retries << '\n' << result.cigar << '\n';
            return 0;
        }
    }
    return 1;
}
