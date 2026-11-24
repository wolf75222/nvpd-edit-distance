// Replace only the two allocation API calls in the included production kernels.
#include <cuda_runtime.h>
#include "chps/baselines.hpp"
#include "chps/fasta.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace probe {
using Clock = std::chrono::steady_clock;
struct Slot {
    int* ptr = nullptr;
    std::size_t bytes = 0;
    bool active = false;
};
struct Batch {
    std::array<Slot, 2> slots{};
    bool reuse = false;
    unsigned allocations = 0, frees = 0, hits = 0;
    std::size_t resident = 0, peak_resident = 0, live = 0, peak_live = 0, allocated_bytes = 0;
    double malloc_seconds = 0, free_seconds = 0;

    cudaError_t release_physical(Slot& s) {
        const auto start = Clock::now();
        auto rc = cudaFree(s.ptr);
        free_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        if (rc == cudaSuccess) {
            resident -= s.bytes;
            ++frees;
            s = {};
        }
        return rc;
    }

    cudaError_t acquire(int** ptr, std::size_t bytes) {
        // The production routine requests its index and score buffers in this order.
        for (auto& s : slots) {
            if (s.active) continue;
            if (s.ptr && s.bytes != bytes) {
                auto rc = release_physical(s);
                if (rc != cudaSuccess) return rc;
            }
            if (s.ptr) {
                ++hits;
            } else {
                const auto start = Clock::now();
                auto rc = cudaMalloc(&s.ptr, bytes);
                malloc_seconds += std::chrono::duration<double>(Clock::now() - start).count();
                if (rc != cudaSuccess) return rc;
                s.bytes = bytes;
                resident += bytes;
                allocated_bytes += bytes;
                ++allocations;
                peak_resident = std::max(peak_resident, resident);
            }
            live += bytes;
            peak_live = std::max(peak_live, live);
            s.active = true;
            *ptr = s.ptr;
            return cudaSuccess;
        }
        return cudaErrorMemoryAllocation;
    }

    cudaError_t give_back(int* ptr) {
        for (auto& s : slots) {
            if (!s.active || s.ptr != ptr) continue;
            s.active = false;
            live -= s.bytes;
            return reuse ? cudaSuccess : release_physical(s);
        }
        return cudaErrorInvalidValue;
    }
    void finish() {
        for (auto& s : slots) {
            if (s.active) throw std::logic_error("kernel returned with a live buffer");
            if (s.ptr && release_physical(s) != cudaSuccess) throw std::runtime_error("cudaFree failed");
        }
        if (resident || live || allocations != frees) throw std::logic_error("unbalanced allocation count");
    }
};
Batch* batch;
cudaError_t acquire(int** ptr, std::size_t bytes) { return batch->acquire(ptr, bytes); }
cudaError_t release(int* ptr) { return batch->give_back(ptr); }
} // namespace probe

#define cudaMalloc ::probe::acquire
#define cudaFree ::probe::release
#include "../src/cuda.cu"
#undef cudaMalloc
#undef cudaFree

int main(int argc, char** argv) try {
    if (argc != 3) return 2;
    const auto sequences = chps::read_fasta(argv[1]);
    if (sequences.size() < 2) throw std::invalid_argument("two or more records required");
    const std::string mode = argv[2];
    if (mode != "normal" && mode != "reuse") throw std::invalid_argument("mode must be normal or reuse");
    std::vector<int> gold, distances;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        for (std::size_t j = i + 1; j < sequences.size(); ++j) {
            gold.push_back(chps::edlib_distance(sequences[i].bases, sequences[j].bases).distance);
        }
    }
    distances.reserve(gold.size());
    probe::Batch warm;
    probe::batch = &warm;
    const auto first = chps::nvpd_cuda(sequences[0].bases, sequences[1].bases);
    warm.finish();
    if (first.distance != gold.front()) throw std::runtime_error("warm-up distance mismatch");

    probe::Batch batch;
    batch.reuse = mode == "reuse";
    probe::batch = &batch;
    std::size_t max_workspace = 0;
    const auto start = probe::Clock::now();
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        for (std::size_t j = i + 1; j < sequences.size(); ++j) {
            const auto r = chps::nvpd_cuda(sequences[i].bases, sequences[j].bases);
            distances.push_back(r.distance);
            max_workspace = std::max(max_workspace, r.device_bytes);
        }
    }
    batch.finish();
    const double seconds = std::chrono::duration<double>(probe::Clock::now() - start).count();
    if (gold != distances || max_workspace != batch.peak_live)
        throw std::runtime_error("distance or explicit workspace mismatch");
    std::uint64_t checksum = 14695981039346656037ull;
    for (int d : distances) checksum = (checksum ^ static_cast<std::uint64_t>(d)) * 1099511628211ull;
    std::cout << std::setprecision(12) << "{\"pairs\":" << gold.size() << ",\"seconds\":" << seconds
              << ",\"malloc_calls\":" << batch.allocations << ",\"free_calls\":" << batch.frees
              << ",\"reuse_hits\":" << batch.hits << ",\"allocated_bytes\":" << batch.allocated_bytes
              << ",\"peak_reserved_bytes\":" << batch.peak_resident << ",\"peak_live_bytes\":" << batch.peak_live
              << ",\"malloc_seconds\":" << batch.malloc_seconds << ",\"free_seconds\":" << batch.free_seconds
              << ",\"distance_fnv64\":\"" << std::hex << checksum << "\"}\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}
