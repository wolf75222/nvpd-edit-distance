#include "cigar_replay.hpp"
#include "corridor_plan_check.hpp"
#include "packed_index.hpp"
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
    check_corridor_plan();
    // Vary the position of the last match across the checkpoint and the
    // byte escape, including symbols that have not yet appeared.
    for (int first : {1, 254, 255, 256, 257, 510, 511, 512}) {
        std::string target(1029, 'A');
        target[first - 1] = 'C';
        target[first + 255] = 'G';
        const detail::MatchIndex index("ACGT", target);
        const auto words = detail::pack_match_offsets(index, static_cast<int>(target.size()));
        const auto checkpoints = detail::match_checkpoints(index, static_cast<int>(target.size()));
        for (int j = 0; j <= static_cast<int>(target.size()); ++j)
            for (int symbol = 0; symbol < 4; ++symbol)
                if (detail::unpack_match(words[j], checkpoints.data(), j, symbol) !=
                    index.entries[symbol * words.size() + j])
                    throw std::runtime_error("packed checkpoint differs from the full match index");
    }

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
    // Absent symbols and a match preceding the whole block must also work on SIMD tails.
    const std::string query = "RNNNAYNNNNNNNNNR";
    for (int length : {16, 17, 31, 32, 33, 63, 64, 65}) {
        std::string target(length, 'A');
        target.front() = 'R';
        target[15] = 'Y';
        target.back() = 'G';
        const int expected = sequential(query, target, false).distance;
        if (nvpd_omp(query, target).distance != expected || nvpd_omp(target, query).distance != expected)
            throw std::runtime_error("wrong NvPD score at a SIMD block boundary");
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
    // Every block has a short but distinct halo. Exercise both full time
    // tiles and an incomplete final tile after the last match in a block.
    for (int rows : {16, 17, 31, 32, 33, 65}) {
        std::string query(rows, 'A'), target(1039, 'A');
        std::mt19937 engine(47);
        for (char& c : query) c = "ACGT"[engine() % 4];
        for (char& c : target) c = "ACGT"[engine() % 4];
        const int reference = edlib_distance(query, target).distance;
        if (nvpd_cuda(query, target).distance != reference || nvpd_cuda(target, query).distance != reference)
            throw std::runtime_error("wrong compact CUDA score at a time-tile boundary");
    }
    {
        std::mt19937 engine(3701);
        std::string query(32801, 'A');
        for (char& c : query) c = "ACGT"[engine() % 4];
        // An insertion followed much later by a deletion moves the optimal
        // path to the edge of a narrow band through many complete time tiles.
        std::string target = query;
        target.insert(31, 1, 'T');
        target.erase(31999, 1);
        target[16512] = target[16512] == 'A' ? 'C' : 'A';
        CudaWorkspace workspace;
        for (int inserted : {0, 1, 127}) {
            target.append(inserted, 'G');
            const int expected = edlib_distance(query, target).distance;
            if (workspace.distance(query, target).distance != expected ||
                workspace.distance(target, query).distance != expected ||
                nvpd_cuda_persistent(query, target).distance != expected ||
                nvpd_cuda_persistent(target, query).distance != expected)
                throw std::runtime_error("wrong CUDA score at an offset corridor edge");
        }
        replay("ACT", "ACCT", nvpd_cuda_persistent("ACT", "ACCT", true));
        if (nvpd_cuda_persistent("", "ACCT").distance != 4 || nvpd_cuda_persistent(query, query).distance != 0)
            throw std::runtime_error("wrong persistent CUDA boundary score");
        query[1919] = 'N';
        if (workspace.distance(query, target).distance != edlib_distance(query, target).distance ||
            nvpd_cuda_persistent(query, target).distance != edlib_distance(query, target).distance)
            throw std::runtime_error("wrong literal ambiguity score after a CUDA corridor call");
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
