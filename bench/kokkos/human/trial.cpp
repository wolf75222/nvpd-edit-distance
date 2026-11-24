#include "chps/dp.hpp"
#include "chps/fasta.hpp"
#include "chps/kokkos.hpp"
#include "row_policy.hpp"
#include <omp.h>
#ifdef CHPS_HUMAN_PROFILE
#include <cuda_profiler_api.h>
#endif
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <sys/resource.h>

using namespace chps;
using View = std::string_view;
void check() {
    std::mt19937 rng(2026);
    for (int id = 0; id < 101; ++id) {
        std::string a(id % 30, 'A'), b((17 * id) % 53, 'N');
        for (char &c : a)
            c = "ACGTNRY"[rng() % 7];
        for (char &c : b)
            c = "ACGTNRY"[rng() % 7];
        if (id == 80) {
            a = std::string(513, 'A') + "C";
            b = "G" + std::string(512, 'A');
        }
        const auto gold = sequential(a, b, true);
        for (auto implementation : {detail::nvpd_kokkos_rows, detail::nvpd_kokkos_opt, nvpd_kokkos}) {
            const auto cost = implementation(a, b, false);
            const auto path = implementation(a, b, true);
            if (cost.distance != gold.distance || path.distance != gold.distance ||
                path.cigar != gold.cigar)
                throw std::runtime_error("cost or traceback did not match the sequential reference");
        }
    }
}

int main(int argc, char **argv) try {
    if (argc == 2 && View(argv[1]) == "--check") {
        KokkosSession session(2);
        check();
        return 0;
    }
    if (argc != 5)
        return 2;
    const View algorithm = argv[1];
    const int workers = std::stoi(argv[3]), expected = std::stoi(argv[4]);
    omp_set_dynamic(0);
    omp_set_num_threads(workers);
    const auto seq = read_fasta(argv[2]);
    if (seq.size() != 2)
        throw std::runtime_error("expected one sequence pair");
    std::optional<KokkosSession> runtime;
    auto kernel = nvpd_omp;
    const char *space = "OpenMP";
    if (algorithm == "native") {
#ifdef CHPS_CUDA
        kernel = nvpd_cuda;
        space = "Cuda";
#endif
    } else {
        runtime.emplace(workers);
        space = kokkos_execution_space();
        if (algorithm == "kokkos-rows")
            kernel = detail::nvpd_kokkos_rows;
        else if (algorithm == "kokkos-opt")
            kernel = detail::nvpd_kokkos_opt;
        else
            throw std::invalid_argument("unknown algorithm");
    }
    if (kernel(seq[0].bases, seq[1].bases, false).distance != expected)
        throw std::runtime_error("warmup mismatch");
#ifdef CHPS_HUMAN_PROFILE
    if (cudaProfilerStart() != cudaSuccess)
        throw std::runtime_error("cannot start CUDA profiler range");
#endif
    const auto start = std::chrono::steady_clock::now();
    const auto result = kernel(seq[0].bases, seq[1].bases, false);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
#ifdef CHPS_HUMAN_PROFILE
    if (cudaProfilerStop() != cudaSuccess)
        throw std::runtime_error("cannot end CUDA profiler range");
#endif
    if (result.distance != expected)
        throw std::runtime_error("timed score mismatch");
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage))
        return 2;
    std::cout << space << ',' << result.distance << ',' << std::setprecision(12) << seconds << ','
              << usage.ru_maxrss << ',' << result.device_bytes << '\n';
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
