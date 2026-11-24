#include "chps/baselines.hpp"
#include "chps/cuda.hpp"
#include "chps/fasta.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cuda_profiler_api.h>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string_view>
#include <sys/resource.h>

extern int roof_variant;
using namespace chps;
using View = std::string_view;
#ifdef CHPS_ROOF_LOCAL
namespace chps {
Result experiment_cuda(View, View, bool);
}
#ifdef CHPS_ROOF_CLASSES
namespace chps {
Result class_local(View, View, bool);
}
#endif

static void verify_local() {
    CudaWorkspace workspace(CudaHalo::local);
    for (int n : {0, 1, 31, 255, 256, 257, 511, 513, 4097, 65011, 301001, 65, 4097}) {
        std::string b(n, 'A');
        for (int j = 0; j < n; j += 7)
            b[j] = "CGT"[j % 3];
        if (n > 1) {
            b[0] = 'R';
            b[n - 1] = 'R';
        }
        for (int m : {0, 1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 24, 33}) {
            std::string a(m, 'A');
            for (int i = 0; i < m; i += 3)
                a[i] = "RCGT"[i % 4];
            const auto expected = edlib_distance(a, b).distance;
            if (workspace.distance(a, b).distance != expected ||
                workspace.distance(b, a).distance != expected)
                throw std::runtime_error("local tile score mismatch");
#ifdef CHPS_ROOF_CLASSES
            if (class_local(a, b, false).distance != expected ||
                class_local(b, a, false).distance != expected)
                throw std::runtime_error("two-class score mismatch");
#endif
            if (n <= 513) {
                auto trace = workspace.distance(a, b, true), gold = sequential(a, b, true);
                if (trace.distance != expected || trace.cigar != gold.cigar)
                    throw std::runtime_error("trace mismatch");
            }
        }
    }
}
#endif

static void verify() {
    std::mt19937 rng(20261124);
    CudaWorkspace work;
    for (int t = 0; t < 102; ++t) {
        const int m = t % 15, n = std::array{0,  1,  2,  15,  16,  31,  32,  33,
                                             63, 64, 65, 255, 256, 257, 511, 513}[t % 16];
        std::string a(m, 'A'), b(n, 'A');
        for (char &c : a)
            c = "ACGTRYN"[rng() % 7];
        for (char &c : b)
            c = "ACGTRYN"[rng() % 7];
        if (n > 64 && t % 3 == 0) {
            // Put the distant predecessor across the warp and block boundaries.
            b.assign(n, 'C');
            b[0] = 'R';
            b[n - 1] = 'R';
            a.assign(m, 'R');
        }
        if (t % 2)
            std::swap(a, b);
        auto gold = sequential(a, b, true);
        for (roof_variant = 0; roof_variant < 3; ++roof_variant) {
#ifdef CHPS_ROOF_LOCAL
            auto cost = experiment_cuda(a, b, false), trace = experiment_cuda(a, b, true);
#else
            auto cost = work.distance(a, b, false), trace = work.distance(a, b, true);
#endif
            if (cost.distance != gold.distance || trace.distance != gold.distance ||
                trace.cigar != gold.cigar)
                throw std::runtime_error("differential CUDA score/trace failure");
        }
    }
    std::cout << "102 pairs x three variants: scalar score and trace matched, including empty/odd and "
                 "warp tails\n";
}

int main(int argc, char **argv) try {
    if (argc == 2 && View(argv[1]) == "--check") {
        verify();
#ifdef CHPS_ROOF_LOCAL
        verify_local();
        std::cout << "195 local-halo sizes in both orientations; empty/tile remainders/long dependencies "
                     "and reuse passed\n";
#endif
        return 0;
    }
    if (argc != 5)
        throw std::invalid_argument("METHOD INPUT.fa EXPECTED profile|time");
    View mode = argv[1];
    auto records = read_fasta(argv[2]);
    if (records.size() != 2)
        throw std::invalid_argument("expected two records");
    const auto &a = records[0].bases;
    const auto &b = records[1].bases;
    const int expected = std::stoi(argv[3]);
    bool profile = View(argv[4]) == "profile";
    if (mode == "native")
        roof_variant = 0;
    else if (mode == "fused")
        roof_variant = 1;
    else if (mode == "shuffle")
        roof_variant = 2;
#ifdef CHPS_ROOF_LOCAL
    else if (mode == "local")
        roof_variant = 0;
#endif
#ifdef CHPS_ROOF_CLASSES
    else if (mode == "classes")
        roof_variant = 0;
#endif
    else if (mode != "edlib" && mode != "wfa2")
        throw std::invalid_argument("unknown method");
    auto score = [&] {
        if (mode == "edlib")
            return edlib_distance(a, b);
        if (mode == "wfa2")
            return wfa_distance(a, b);
#ifdef CHPS_ROOF_LOCAL
        if (mode == "local")
            return nvpd_cuda_local(a, b);
        if (mode == "fused" || mode == "shuffle")
            return experiment_cuda(a, b, false);
#endif
#ifdef CHPS_ROOF_CLASSES
        if (mode == "classes")
            return class_local(a, b, false);
#endif
        return nvpd_cuda(a, b, false);
    };
    auto warm = score();
    if (warm.distance != expected)
        throw std::runtime_error("warm score mismatch");
    if (profile && cudaProfilerStart() != cudaSuccess)
        throw std::runtime_error("cudaProfilerStart failed");
    const auto tick = std::chrono::steady_clock::now();
    auto timed = score();
    const auto dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - tick).count();
    if (profile && cudaProfilerStop() != cudaSuccess)
        throw std::runtime_error("cudaProfilerStop failed");
    if (timed.distance != expected)
        throw std::runtime_error("timed score mismatch");
    rusage rss{};
    if (getrusage(RUSAGE_SELF, &rss))
        throw std::runtime_error("getrusage failed");
    std::cout << a.size() << ',' << b.size() << ',' << timed.distance << ',' << std::setprecision(12)
              << dt << ',' << rss.ru_maxrss << ',' << timed.device_bytes << '\n';
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
