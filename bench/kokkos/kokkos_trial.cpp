// One process per trial so the host high-water mark belongs to one backend.
#include "chps/fasta.hpp"
#include "chps/kokkos.hpp"
#include <omp.h>
#include <sys/resource.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    try {
        const std::string algorithm = argv[1];
        const int threads = std::stoi(argv[3]);
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
        auto kernel = chps::nvpd_omp;
        const char* space = "OpenMP";
        std::optional<chps::KokkosSession> runtime;
        if (algorithm == "nvpd-kokkos") {
            runtime.emplace(threads);
            kernel = chps::nvpd_kokkos;
            space = chps::kokkos_execution_space();
        } else if (algorithm == "nvpd-cuda") {
#ifdef CHPS_CUDA
            kernel = chps::nvpd_cuda;
            space = "CUDA";
#else
            throw std::runtime_error("CUDA is disabled");
#endif
        } else if (algorithm != "nvpd-omp") {
            throw std::invalid_argument("unknown benchmark backend");
        }
        const auto records = chps::read_fasta(argv[2]);
        if (records.size() != 2) throw std::invalid_argument("expected two FASTA records");
        const auto& a = records[0].bases;
        const auto& b = records[1].bases;
        const int check = kernel(a, b, false).distance;
        const auto start = std::chrono::steady_clock::now();
        const auto res = kernel(a, b, false);
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        if (res.distance != check) throw std::runtime_error("warm and measured distances differ");
        rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage) != 0) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
        usage.ru_maxrss /= 1024;
#endif
        std::cout << std::setprecision(12) << "{\"execution_space\":\"" << space
                  << "\",\"m\":" << a.size() << ",\"n\":" << b.size()
                  << ",\"distance\":" << res.distance << ",\"seconds\":" << elapsed.count()
                  << ",\"peak_rss_kib\":" << usage.ru_maxrss << ",\"device_bytes\":" << res.device_bytes << "}\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }
}
