#include "cigar_replay.hpp"
#include "chps/baselines.hpp"
#include "chps/banded.hpp"
#ifdef CHPS_CUDA
#include "chps/cuda.hpp"
#endif
#ifdef CHPS_KOKKOS
#include "chps/kokkos.hpp"
#endif
#include <omp.h>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace chps;

static void check(const std::string& a, const std::string& b) {
    int reference = sequential(a, b, false).distance;
    std::vector<Result> costs{antidiagonal(a,b), nvpd_omp(a,b), edlib_distance(a,b), wfa_distance(a,b), biwfa(a,b)};
    costs.push_back(banded(a,b));
    costs.push_back(banded_trim(a,b));
    std::vector<Result> traces{edlib_trace(a,b), wfa_trace(a,b), biwfa(a,b,true), hirschberg(a,b,true), sequential(a,b,true), nvpd_omp(a,b,true)};
#ifdef CHPS_KOKKOS
    costs.push_back(nvpd_kokkos(a,b));
    traces.push_back(nvpd_kokkos(a,b,true));
#endif
#ifdef CHPS_CUDA
    costs.push_back(nvpd_cuda(a,b));
    traces.push_back(nvpd_cuda(a,b,true));
    costs.push_back(nvpd_cuda_local(a,b));
    traces.push_back(nvpd_cuda_local(a,b,true));
#endif
    for (const auto& r : traces) replay(a, b, r);
    costs.insert(costs.end(), traces.begin(), traces.end());
    for (const auto& r : costs)
        if (r.distance != reference) throw std::runtime_error("wrong distance: " + a + " / " + b);
}

int main() {
    omp_set_num_threads(3);
#ifdef CHPS_KOKKOS
    KokkosSession session(3);
#endif
    std::vector<std::string> words{"", "A", "AA", "G", "GA", "N", "AN", "AGT", "TTTT", "ACGTN", "RYNN"};
    for (const auto& a : words) for (const auto& b : words) check(a,b);
    std::mt19937 rng(2026);
    auto dna = [&](unsigned n) {
        std::string s;
        while (s.size() < n) s += "ACGTNR"[rng() % 6];
        return s;
    };
    for (int i = 0; i < 70; ++i) {
        auto a = dna(rng() % 65);
        auto b = dna(rng() % 90);
        check(a,b);
    }
    check(dna(259), dna(519));
    // Force bidirectional breakpoints beyond upstream's 250-edit base cases.
    const std::string left = dna(917), right = dna(913);
    for (const auto& path : {biwfa(left, right, true), wfa_trace(left, right)}) {
        replay(left, right, path);
        if (path.distance != edlib_distance(left, right).distance ||
            path.distance != biwfa(left, right).distance || path.distance <= 250)
            throw std::runtime_error("wrong wavefront breakpoint or score");
    }
    // Exercise absent symbols, distant matches and SIMD tails in both orientations.
    const std::string rare = "R" + std::string(511, 'A') + "Y";
    check("RYNNR", rare);
    check(rare, "RYNNR");
    check("ARNNNRRYYYY", rare);
    std::string repeated;
    for (int copy = 0; copy < 85; ++copy) repeated += "ACGT";
    check("ACNTGACGTNACG", repeated);
    for (bool transpose : {false, true}) {
        std::string a = dna(20'011), b = dna(37);
        if (transpose) std::swap(a, b);
        const auto path = hirschberg(a, b, true);
        replay(a, b, path);
        if (path.distance != edlib_distance(a, b).distance)
            throw std::runtime_error("wrong linear-space traceback beyond dense limit");
    }
    // One DP row isolates column indexing, including million-column inputs.
    for (int length : {127, 128, 129, 254, 255, 256, 257, 65535, 1000000}) {
        const std::string absent(length, 'A');
        if (nvpd_omp("R", absent).distance != length)
            throw std::runtime_error("wrong NvPD distance for absent symbol");
    }
#ifdef CHPS_CUDA
    {
        CudaWorkspace workspace(CudaHalo::local);
        std::string target(5179, 'A');
        for (std::size_t j = 0; j + 3 <= target.size(); j += 263) target.replace(j, 3, "CGT");
        const std::string query = "CGTACGTACGTACGTAC";
        if (workspace.distance(query, target).distance != edlib_distance(query, target).distance)
            throw std::runtime_error("wrong local CUDA score across packed-offset escape and tile tails");
    }
    {
        CudaWorkspace workspace;
        const std::string first(530, 'A'), second(500, 'A');
        const auto path = workspace.distance(std::string(17, 'A'), first, true);
        replay(std::string(17, 'A'), first, path);
        const auto cost = workspace.distance("ACGT", second);
        const auto expected_peak = (18 * (first.size() + 1) + 4 * (second.size() + 1)) * sizeof(int);
        const auto idle = workspace.distance("", second);
        if (cost.distance != edlib_distance("ACGT", second).distance ||
            cost.device_bytes != expected_peak || idle.device_bytes != 6 * (second.size() + 1) * sizeof(int))
            throw std::runtime_error("wrong score or allocation accounting after CUDA workspace resize");
    }
    for (int spacing : {4, 64, 65}) {
        std::string target(400003, 'A');
        for (std::size_t j = 0; j < target.size(); j += spacing) target[j] = 'C';
        std::string query = "AACAAACAAAAACACAA";
        if (nvpd_cuda(query, target).distance != edlib_distance(query, target).distance)
            throw std::runtime_error("wrong CUDA distance at wide block or halo boundary");
    }
#endif
    if (sequential("GATTACA", "GCATGCU").distance != 4) return 1;
    std::cout << "196 differential pairs with traceback replay and 9 NvPD column widths\n";
}
