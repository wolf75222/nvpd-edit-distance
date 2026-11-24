#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/split.hpp"
#include "cigar_replay.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
using Method = chps::Result (*)(std::string_view, std::string_view, bool, chps::SplitStats&, std::size_t);

chps::Result front(std::string_view a, std::string_view b, bool trace, chps::SplitStats& stats, std::size_t) {
    ++stats.wavefront_pieces;
    return chps::biwfa(a, b, trace);
}

template<chps::SplitPolicy policy>
chps::Result split(std::string_view a, std::string_view b, bool trace, chps::SplitStats& stats, std::size_t lower) {
    return chps::nvpd_split_profile(a, b, trace, stats, policy, lower);
}

struct Cell {
    chps::Result result{};
    chps::SplitStats stats;
};

void measure(const std::vector<chps::Record>& records, const std::string& oracle,
             std::string_view name, Method solve, bool trace, int threads, int block, int position) {
    const int n = static_cast<int>(records.size()), groups = std::min(threads, n - 1);
    std::vector<Cell> matrix(n * n);
    std::exception_ptr failure;
    const auto start = std::chrono::steady_clock::now();
    std::vector<chps::SplitComposition> composition;
    if (name == "split-filter") {
        composition.reserve(n);
        for (const auto& record : records) composition.emplace_back(record.bases);
    }
    #pragma omp parallel num_threads(threads)
    {
        omp_set_num_threads(1);
        #pragma omp for collapse(2) schedule(dynamic, 1)
        for (int i = 0; i < n; ++i) {
            for (int group = 0; group < groups; ++group) {
                try {
                    for (int j = i + 1 + group; j < n; j += groups) {
                        auto& cell = matrix[i*n+j];
                        const auto lower = composition.empty() ? 0 : composition[i].lower_bound(composition[j]);
                        cell.result = solve(records[i].bases, records[j].bases, trace, cell.stats, lower);
                    }
                } catch (...) {
                    #pragma omp critical
                    { if (!failure) failure = std::current_exception(); }
                }
            }
        }
    }
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    if (failure) std::rethrow_exception(failure);
    std::ostringstream tsv;
    tsv << "id";
    for (const auto& record : records) tsv << '\t' << record.id;
    tsv << '\n';
    int proposed = 0, accepted = 0, with_proposals = 0, with_cuts = 0, retries = 0;
    int dense = 0, fronts = 0, with_dense = 0, with_fronts = 0, pruned = 0;
    for (int i = 0; i < n; ++i) {
        tsv << records[i].id;
        for (int j = 0; j < n; ++j) {
            tsv << '\t' << matrix[std::min(i,j)*n+std::max(i,j)].result.distance;
            if (i >= j) continue;
            if (trace) replay(records[i].bases, records[j].bases, matrix[i*n+j].result);
            const auto& s = matrix[i*n+j].stats;
            proposed += s.proposed_cuts;
            accepted += s.accepted_cuts;
            with_proposals += s.proposed_cuts > 0;
            with_cuts += s.accepted_cuts > 0;
            retries += s.global_retry;
            pruned += s.composition_pruned;
            dense += s.dense_pieces;
            fronts += s.wavefront_pieces;
            with_dense += s.dense_pieces > 0;
            with_fronts += s.wavefront_pieces > 0;
        }
        tsv << '\n';
    }
    if (tsv.str() != oracle) throw std::runtime_error("distance matrix differs from reference");
    std::cout << block << ',' << position << ',' << name << ',' << (trace ? "trace" : "score")
              << ',' << threads << ',' << n*(n-1)/2
              << ',' << std::setprecision(10) << elapsed.count() << ',' << proposed << ',' << accepted
              << ',' << with_proposals << ',' << with_cuts << ',' << retries << ',' << dense << ',' << fronts
              << ',' << with_dense << ',' << with_fronts << ',' << pruned << '\n';
}
}

int main(int argc, char** argv) {
    if (argc != 4) throw std::invalid_argument("ablation COHORT.fa MATRIX.tsv THREADS");
    const auto records = chps::read_fasta(argv[1]);
    const int threads = std::stoi(argv[3]);
    if (records.size() != 64 || threads < 1) throw std::invalid_argument("expected 64 sequences and positive threads");
    std::ifstream ref(argv[2], std::ios::binary);
    if (!ref) throw std::runtime_error("cannot open reference matrix");
    const std::string oracle{std::istreambuf_iterator<char>(ref), {}};
    constexpr std::array<std::string_view,4> names{"biwfa", "trim-select", "split", "split-filter"};
    constexpr std::array<Method,4> methods{front, split<chps::SplitPolicy::single_piece>,
               split<chps::SplitPolicy::retain_verified>, split<chps::SplitPolicy::retain_verified>};
    // Each method occupies all four positions in measured blocks 1–4.
    constexpr int orders[][4]{{0,1,2,3}, {0,1,3,2}, {1,2,0,3}, {2,3,1,0}, {3,0,2,1}, {2,0,3,1}};
    std::cout << "block,position,method,mode,threads,pairs,seconds,proposed_cuts,accepted_cuts,"
                 "pairs_proposed,pairs_accepted,pairs_global_retry,dense_calls,front_calls,"
                 "pairs_dense,pairs_front,pairs_composition_pruned\n";
    for (bool trace : {false,true}) {
        int block = 0;
        for (const auto& order : orders) {
            int position = 0;
            for (int m : order) measure(records, oracle, names[m], methods[m], trace, threads, block, ++position);
            ++block;
        }
    }
}
