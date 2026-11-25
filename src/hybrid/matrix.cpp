#include "matrix.hpp"
#include "sketch.hpp"
#include "chps/baselines.hpp"
#ifdef CHPS_CUDA
#include "chps/cuda.hpp"
#endif

#include <atomic>
#include <chrono>
#include <exception>
#include <mutex>
#include <omp.h>
#include <stdexcept>

namespace chps::hybrid {
namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

enum class Solver { edlib, wfa2, cuda };
struct Pair {
    std::size_t first, second;
    Solver solver;
};

Solver choose(Policy policy, const std::vector<Sketch>& words, std::size_t i,
              std::size_t j, bool even) {
    switch (policy) {
    case Policy::edlib: return Solver::edlib;
    case Policy::wfa2: return Solver::wfa2;
    case Policy::gpu: return Solver::cuda;
    case Policy::static_split: return even ? Solver::cuda : Solver::edlib;
    case Policy::hybrid:
    case Policy::hybrid_share:
    case Policy::cpu_sketch:
        if (words[i].empty() || words[j].empty()) return Solver::edlib;
        if (overlap(words[i], words[j]).close()) return Solver::wfa2;
        return policy == Policy::cpu_sketch ? Solver::edlib : Solver::cuda;
    }
    throw std::logic_error("invalid routing policy");
}
}

Matrix align(const std::vector<Record>& records, Policy policy, int threads) {
    if (threads <= 0) throw std::invalid_argument("threads must be positive");
#ifndef CHPS_CUDA
    if (policy != Policy::edlib && policy != Policy::wfa2 && policy != Policy::cpu_sketch)
        throw std::invalid_argument("this policy requires a build with CHPS_CUDA=ON");
#endif
    auto start = Clock::now();
    Matrix result;
    const auto n = records.size();
    if (n && n > result.scores.max_size() / n) throw std::length_error("matrix too large");
    result.scores.resize(n * n);
    std::vector<Pair> cpu, gpu;
    std::exception_ptr failure;
    std::vector<Sketch> words;
    const bool share = policy == Policy::hybrid_share;
    const bool estimate = share || policy == Policy::hybrid || policy == Policy::cpu_sketch;
    if (estimate) {
        const auto sketch_start = Clock::now();
        words.resize(n);
        #pragma omp parallel for num_threads(threads) schedule(static)
        for (std::size_t i = 0; i < n; ++i) {
            try {
                words[i] = sketch(records[i].bases);
            } catch (...) {
                #pragma omp critical(chps_route_error)
                { if (!failure) failure = std::current_exception(); }
            }
        }
        if (failure) std::rethrow_exception(failure);
        result.sketch_seconds = elapsed(sketch_start);
    }
    const auto route_start = Clock::now();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const auto solver = choose(policy, words, i, j, (gpu.size() + cpu.size()) % 2 == 0);
            (solver == Solver::cuda ? gpu : cpu).push_back({i, j, solver});
        }
    }
    if (estimate) result.sketch_seconds += elapsed(route_start);
    std::size_t front = 0, back = gpu.size();
    std::mutex queue_lock;
    auto take_far = [&](bool from_front) -> const Pair* {
        const std::lock_guard lock(queue_lock);
        if (front == back) return nullptr;
        return &gpu[from_front ? front++ : --back];
    };
    std::atomic<std::size_t> next{0};
    auto store = [&](const Pair& p, int score) {
        result.scores[p.first * n + p.second] = score;
        result.scores[p.second * n + p.first] = score;
    };
    const int workers = cpu.empty() && !share ? 1 : threads;
    #pragma omp parallel num_threads(workers)
    {
        try {
#ifdef CHPS_CUDA
            if (omp_get_thread_num() == 0 && !gpu.empty()) {
                CudaWorkspace workspace(CudaHalo::local);
                auto run_gpu = [&](const Pair& p) {
                    store(p, workspace.distance(records[p.first].bases, records[p.second].bases).distance);
                    ++result.gpu_pairs;
                };
                if (share) {
                    while (const auto* p = take_far(true)) run_gpu(*p);
                } else {
                    for (const auto& p : gpu) run_gpu(p);
                }
            }
#endif
            if (gpu.empty() || omp_get_thread_num() != 0 || omp_get_num_threads() == 1) {
                for (;;) {
                    auto index = next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= cpu.size()) break;
                    const auto& p = cpu[index];
                    const auto kernel = p.solver == Solver::edlib ? edlib_distance : wfa_distance;
                    store(p, kernel(records[p.first].bases, records[p.second].bases).distance);
                }
                if (share) {
                    while (const auto* p = take_far(false))
                        store(*p, edlib_distance(records[p->first].bases, records[p->second].bases).distance);
                }
            }
        } catch (...) {
            #pragma omp critical(chps_route_error)
            { if (!failure) failure = std::current_exception(); }
        }
    }
    if (failure) std::rethrow_exception(failure);
    result.cpu_pairs = cpu.size() + gpu.size() - result.gpu_pairs;
    result.seconds = elapsed(start);
    return result;
}

} // namespace chps::hybrid
