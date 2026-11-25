#include "chps/fasta.hpp"
#include "chps/split.hpp"

#include <array>
#include <chrono>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
constexpr int workers = 4, repetitions = 16;

struct Pair {
    std::string path;
    std::size_t m{}, n{};
    int score{};
    std::vector<chps::Record> records;
};

void load(Pair& p) {
    p.records = chps::read_fasta(p.path);
    if (p.records.size() != 2 || p.records[0].bases.size() != p.m || p.records[1].bases.size() != p.n)
        throw std::runtime_error("FASTA lengths differ from reference: " + p.path);
}

template<class F> double parallel(F work, std::array<int, workers>& cpus) {
    std::array<std::exception_ptr, workers> errors;
    const auto start = Clock::now();
    #pragma omp parallel num_threads(workers)
    {
        const int tid = omp_get_thread_num();
        try {
            omp_set_num_threads(1);
            const int place = omp_get_place_num();
            if (omp_get_num_threads() != workers || place < 0 || omp_get_place_num_procs(place) != 1)
                throw std::runtime_error("Expected four bound workers");
            omp_get_place_proc_ids(place, &cpus[tid]);
            work(tid);
        } catch (...) {
            errors[tid] = std::current_exception();
        }
    }
    const auto stop = Clock::now();
    for (const auto& error : errors)
        if (error) std::rethrow_exception(error);
    return std::chrono::duration<double>(stop - start).count();
}
} // namespace

int main(int argc, char** argv) try {
    if (argc != 3) throw std::invalid_argument("Usage: grace_placement inputs.tsv master|owner");
    const std::string mode = argv[2];
    if (mode != "master" && mode != "owner") throw std::invalid_argument("Unknown placement mode");
    std::ifstream in(argv[1]);
    std::array<Pair, workers> pairs;
    std::string extra;
    for (auto& p : pairs)
        if (!(in >> p.path >> p.m >> p.n >> p.score)) throw std::runtime_error("Expected four chromosomes");
    if (in >> extra) throw std::runtime_error("Surplus chromosomes");
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    omp_set_num_threads(1);
    std::array<int, workers> cpus{};
    const auto start = Clock::now();
    if (mode == "master") {
        for (auto& pair : pairs) load(pair);
    } else {
        parallel([&](int tid) { load(pairs[tid]); }, cpus);
    }
    const double load_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    const auto align = [&](int tid, int count) {
        const auto& p = pairs[tid];
        for (int k = 0; k < count; ++k)
            if (chps::nvpd_sparse(p.records[0].bases, p.records[1].bases).distance != p.score)
                throw std::runtime_error("Score differs from reference: " + p.path);
    };
    const auto warm = parallel([&](int tid) { align(tid, 1); }, cpus);
    std::cout << std::setprecision(10) << "ready " << load_seconds << ' ' << warm;
    for (auto cpu : cpus) std::cout << ' ' << cpu;
    for (const auto& p : pairs)
        for (const auto& rec : p.records)
            std::cout << ' ' << reinterpret_cast<std::uintptr_t>(rec.bases.data()) << ' ' << rec.bases.size();
    std::cout << std::endl;
    std::string go;
    if (!std::getline(std::cin, go) || go != "run") throw std::runtime_error("Missing measurement command");
    std::cout << parallel([&](int tid) { align(tid, repetitions); }, cpus) << std::endl;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
