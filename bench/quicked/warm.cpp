#include "chps/fasta.hpp"
#if defined(CHPS_BENCH_NATIVE)
#include "chps/baselines.hpp"
#elif defined(CHPS_BENCH_ASTAR)
#include "astarpa.h"
#else
#include "chps/quicked.h"
#endif
#include <chrono>
#include <climits>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const auto pair = chps::read_fasta(argv[1]);
    if (pair.size() != 2) return 2;
    const std::string& a = pair[0].bases;
    const std::string& b = pair[1].bases;
    if (a.empty() || b.empty() || std::max(a.size(), b.size()) > INT_MAX / 100 ||
        a.find_first_not_of("ACGT") != a.npos || b.find_first_not_of("ACGT") != b.npos) return 2;
    const std::string mode = argv[2];
    const int calls = std::stoi(argv[3]);
    for (int i = 0; i <= calls; ++i) {
        int score = -1, retries = 0;
        const char* source = mode.c_str();
        const auto start = std::chrono::steady_clock::now();
#if defined(CHPS_BENCH_NATIVE)
        if (mode == "edlib") score = chps::edlib_distance(a,b).distance;
        else if (mode == "wfa2") score = chps::wfa_distance(a,b).distance;
        else return 2;
#elif defined(CHPS_BENCH_ASTAR)
        if (mode != "astarpa2-simple" && mode != "astarpa2-full") return 2;
        auto align = mode == "astarpa2-simple" ? astarpa2_simple : astarpa2_full;
        uint8_t* cigar = nullptr;
        uintptr_t len = 0;
        score = align(reinterpret_cast<const uint8_t*>(a.data()), a.size(),
                      reinterpret_cast<const uint8_t*>(b.data()), b.size(), &cigar, &len);
        if (!cigar) return 1;
        astarpa_free_cigar(cigar);
#elif defined(CHPS_BENCH_STOCK)
        if (mode != "quicked-stock") return 2;
        score = chps_quicked_default(a.data(), a.size(), b.data(), b.size());
#else
        if (mode != "quicked-certified") return 2;
        const auto res = chps_quicked_certified(a.data(), a.size(), b.data(), b.size());
        score = res.distance;
        source = res.source;
        retries = res.retries;
#endif
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        if (i) std::cout << std::setprecision(14) << i << ',' << score << ',' << elapsed.count()
                         << ',' << source << ',' << retries << '\n';
    }
}
