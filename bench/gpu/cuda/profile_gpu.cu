#include "chps/dp.hpp"
#include "chps/fasta.hpp"
#include <cuda_profiler_api.h>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Record the entire backend call. Nsight capture brackets only the last warm call.
int main(int argc, char** argv) try {
    if (argc != 8) throw std::invalid_argument("variant FASTA repeats run trace profile expected");
    const auto records = chps::read_fasta(argv[2]);
    if (records.size() != 2) throw std::invalid_argument("expected a two-record FASTA");
    const auto& a = records[0].bases;
    const auto& b = records[1].bases;
    const int repeats = std::stoi(argv[3]);
    const int run = std::stoi(argv[4]);
    const bool trace = std::stoi(argv[5]) != 0;
    const bool profile = std::stoi(argv[6]) != 0;
    const int expected = std::stoi(argv[7]);
    for (int trial = 0; trial <= repeats; ++trial) {
        const bool capture = profile && trial == repeats;
        if (capture && cudaProfilerStart() != cudaSuccess) throw std::runtime_error("profiler start");
        const auto start = std::chrono::steady_clock::now();
        const auto result = chps::nvpd_cuda(a, b, trace);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (capture && cudaProfilerStop() != cudaSuccess) throw std::runtime_error("profiler stop");
        if (result.distance != expected) throw std::runtime_error("distance differs from validated reference");
        std::cout << argv[1] << "," << std::filesystem::path(argv[2]).filename().string()
                  << "," << a.size() << "," << b.size() << "," << (trial ? "warm" : "cold")
                  << "," << run << "," << trial << "," << trace << "," << std::setprecision(10)
                  << seconds << "," << result.distance << "," << result.device_bytes << "\n";
    }
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return 1;
}
