#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/split.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
constexpr int workers_max = 4;
constexpr int passes_per_worker = 4;

struct Pair {
    std::string a, b;
    int score;
};

std::vector<Pair> load(const char* manifest) {
    std::ifstream in(manifest);
    std::string path;
    std::size_t m, n;
    int score;
    std::vector<Pair> pairs;
    while (in >> path >> m >> n >> score) {
        auto records = chps::read_fasta(path);
        if (records.size() != 2 || records[0].bases.size() != m || records[1].bases.size() != n)
            throw std::runtime_error("Unexpected FASTA records: " + path);
        if (chps::wfa_distance(records[0].bases, records[1].bases).distance != score)
            throw std::runtime_error("Unexpected WFA2 distance: " + path);
        pairs.push_back({std::move(records[0].bases), std::move(records[1].bases), score});
    }
    if (!in.eof() || pairs.size() != 4)
        throw std::runtime_error("Expected four complete chromosome pairs");
    return pairs;
}

struct Measurement {
    double seconds;
    std::array<int, workers_max> cpus;
};

Measurement measure(const std::vector<Pair>& pairs, int workers, int passes) {
    std::array<std::exception_ptr, workers_max> errors{};
    std::array<int, workers_max> cpus{};
    const auto start = Clock::now();
    #pragma omp parallel num_threads(workers)
    {
        const int tid = omp_get_thread_num();
        try {
            omp_set_num_threads(1);
            const int place = omp_get_place_num();
            if (omp_get_num_threads() != workers || place < 0 || omp_get_place_num_procs(place) != 1)
                throw std::runtime_error("Expected one bound hardware thread per OpenMP worker");
            omp_get_place_proc_ids(place, &cpus[tid]);
            for (int pass = 0; pass < passes; ++pass) {
                for (int k = 0; k < static_cast<int>(pairs.size()); ++k) {
                    const auto& pair = pairs[(tid + pass + k) % pairs.size()];
                    if (chps::nvpd_sparse(pair.a, pair.b).distance != pair.score)
                        throw std::runtime_error("NvPD score differs from the chromosome oracle");
                }
            }
        } catch (...) {
            errors[tid] = std::current_exception();
        }
    }
    const auto stop = Clock::now();
    for (const auto& error : errors)
        if (error) std::rethrow_exception(error);
    if (std::set<int>(cpus.begin(), cpus.begin() + workers).size() != static_cast<std::size_t>(workers))
        throw std::runtime_error("Workers share an OpenMP place");
    return {std::chrono::duration<double>(stop - start).count(), cpus};
}
} // namespace

int main(int argc, char** argv) try {
    if (argc != 4) throw std::invalid_argument("Usage: grace_scaling inputs.tsv block permutation");
    const int block = std::stoi(argv[2]);
    const int permutation = std::stoi(argv[3]);
    if (block < 1 || block > 6 || permutation < 0 || permutation > 5)
        throw std::invalid_argument("Invalid block or permutation");
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    omp_set_num_threads(1);
    const auto pairs = load(argv[1]);
    const auto bases = std::accumulate(pairs.begin(), pairs.end(), std::uint64_t{0},
        [](auto sum, const Pair& pair) { return sum + pair.a.size() + pair.b.size(); });
    std::cout << "phase,block,position,mode,workers,pairs,bases,seconds,cpu_ids\n" << std::setprecision(9);
    const auto run = [&](bool warmup, int position, bool strong, int workers, int passes) {
        const auto [seconds, cpus] = measure(pairs, workers, passes);
        const int panel_calls = passes * workers;
        std::cout << (warmup ? "warmup" : "measure") << ',' << block << ',' << position
                  << ',' << (strong ? "strong" : "weak") << ',' << workers
                  << ',' << panel_calls * pairs.size() << ',' << bases * panel_calls << ',' << seconds << ',';
        for (int i = 0; i < workers; ++i) std::cout << (i ? "|" : "") << cpus[i];
        std::cout << std::endl;
    };
    std::array<int, 3> order{1, 2, 4};
    for (int k = 0; k < permutation; ++k) std::next_permutation(order.begin(), order.end());
    for (int i = 0; i < 3; ++i) run(true, i + 1, false, order[i], 1);
    int position = 0;
    for (int workers : order) {
        for (int trial = 0; trial < 2; ++trial) {
            const bool strong = ((block + trial) % 2 == 0);
            const int passes = passes_per_worker * (strong ? workers_max / workers : 1);
            run(false, ++position, strong, workers, passes);
        }
    }
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
