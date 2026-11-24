#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// VmHWM is recorded after exec and includes data, runtime, warmup and timed call.
static long peak_rss_kib() {
    std::ifstream status("/proc/self/status");
    std::string key;
    while (status >> key) {
        if (key == "VmHWM:") {
            long kib;
            if (status >> kib && kib > 0) return kib;
        }
        std::getline(status, key);
    }
    throw std::runtime_error("Linux VmHWM unavailable");
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const auto records = chps::read_fasta(argv[1]);
    const std::string_view algorithm = argv[2], scope = argv[3];
    if (records.size() != 2 || (scope != "score" && scope != "trace")) return 2;
    const auto& a = records[0].bases;
    const auto& b = records[1].bases;
    const bool trace = scope == "trace";
    auto run = [&]() {
        if (algorithm == "biwfa") return chps::biwfa(a, b, trace);
        if (algorithm == "wfa2") return trace ? chps::wfa_trace(a,b) : chps::wfa_distance(a,b);
        if (algorithm == "edlib") return trace ? chps::edlib_trace(a,b) : chps::edlib_distance(a,b);
        throw std::invalid_argument("unknown algorithm");
    };
    const int warm = run().distance;
    const auto start = std::chrono::steady_clock::now();
    const auto result = run();
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    if (result.distance != warm) return 1;
    std::cout << std::setprecision(12) << elapsed.count() << ' ' << peak_rss_kib()
              << ' ' << result.distance << '\n' << result.cigar << '\n';
}
