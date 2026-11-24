#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "cigar_replay.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <utility>

int main(int argc, char** argv) {
    if (argc != 3) throw std::invalid_argument("wfa-trace-batch cohort.fa score|trace");
    const std::string_view mode = argv[2];
    if (mode != "trace" && mode != "score") throw std::invalid_argument("unknown output mode");
    const bool trace = mode == "trace";
    const auto records = chps::read_fasta(argv[1]);
    std::vector<std::pair<int,int>> pairs;
    const int n = static_cast<int>(records.size());
    for (int i = 0; i < n; ++i)
        for (int j = i+1; j < n; ++j) pairs.emplace_back(i,j);
    std::vector<chps::Result> results(pairs.size());
    const auto start = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(dynamic,1)
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        const auto [i,j] = pairs[k];
        results[k] = chps::biwfa(records[i].bases, records[j].bases, trace);
    }
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    std::vector<int> scores(n*n);
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        const auto [i,j] = pairs[k];
        if (trace) replay(records[i].bases, records[j].bases, results[k]);
        scores[i*n+j] = scores[j*n+i] = results[k].distance;
    }
    std::cout << "id";
    for (const auto& record : records) std::cout << '\t' << record.id;
    std::cout << '\n';
    for (int i=0; i<n; ++i) {
        std::cout << records[i].id;
        for (int j=0; j<n; ++j) std::cout << '\t' << scores[i*n+j];
        std::cout << '\n';
    }
    std::cerr << std::setprecision(12) << elapsed.count() << '\n';
}
