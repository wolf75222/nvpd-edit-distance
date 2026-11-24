#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/split.hpp"
#include "cigar_replay.hpp"
#include <sys/resource.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 5) throw std::invalid_argument("trial FASTA algorithm score|trace threads");
    auto records = chps::read_fasta(argv[1]);
    if (records.size() != 2) throw std::invalid_argument("two FASTA records required");
    const auto& a = records[0].bases;
    const auto& b = records[1].bases;
    const std::string algorithm = argv[2];
    const bool trace = std::string_view(argv[3]) == "trace";
    omp_set_num_threads(std::stoi(argv[4]));
    chps::SplitStats stats;
    auto run = [&] {
        if (algorithm == "split-all") return chps::nvpd_split_profile(a,b,trace,stats,chps::SplitPolicy::require_all);
        if (algorithm == "split") return chps::nvpd_split_profile(a,b,trace,stats);
        if (algorithm == "biwfa") return chps::biwfa(a,b,trace);
        if (algorithm == "wfa2") return trace ? chps::wfa_trace(a,b) : chps::wfa_distance(a,b);
        if (algorithm == "edlib") return trace ? chps::edlib_trace(a,b) : chps::edlib_distance(a,b);
        throw std::invalid_argument("unknown method");
    };
    const int warm = run().distance;
    const auto start = std::chrono::steady_clock::now();
    const auto result = run();
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if (warm != result.distance) throw std::runtime_error("warm and measured scores differ");
    if (trace) replay(a,b,result);
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage");
    std::cout << a.size() << ',' << b.size() << ',' << result.distance << ',' << std::setprecision(12)
              << seconds << ',' << usage.ru_maxrss << ',' << stats.proposed_cuts << ',' << stats.accepted_cuts
              << ',' << stats.dense_pieces << ',' << stats.wavefront_pieces << ',' << stats.global_retry << '\n';
}
