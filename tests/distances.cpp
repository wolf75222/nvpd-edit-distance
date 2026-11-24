#include "chps/baselines.hpp"
#include "chps/banded.hpp"
#ifdef CHPS_KOKKOS
#include "chps/kokkos.hpp"
#endif
#include <omp.h>
#include <charconv>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace chps;

static void replay(std::string_view a, std::string_view b, const Result& result) {
    std::string_view cigar = result.cigar;
    size_t i = 0, j = 0;
    int cost = 0;
    while (!cigar.empty()) {
        int count = 0;
        const auto [op, ec] = std::from_chars(cigar.data(), cigar.data() + cigar.size(), count);
        if (ec != std::errc{} || op == cigar.data() + cigar.size() || count <= 0)
            throw std::runtime_error("malformed CIGAR");
        cigar.remove_prefix(op + 1 - cigar.data());
        for (int k = 0; k < count; ++k) {
            switch (*op) {
            case '=': if (a.at(i++) != b.at(j++)) throw std::runtime_error("bad match"); break;
            case 'X': if (a.at(i++) == b.at(j++)) throw std::runtime_error("bad mismatch"); ++cost; break;
            case 'I': (void)b.at(j++); ++cost; break;
            case 'D': (void)a.at(i++); ++cost; break;
            default: throw std::runtime_error("bad CIGAR operation");
            }
        }
    }
    if (i != a.size() || j != b.size() || cost != result.distance)
        throw std::runtime_error("traceback failed to consume sequences or achieve optimal cost");
}

static void check(const std::string& a, const std::string& b) {
    int reference = sequential(a, b, false).distance;
    std::vector<Result> costs{antidiagonal(a,b), nvpd_omp(a,b), edlib_distance(a,b), wfa_distance(a,b)};
    costs.push_back(banded(a,b));
    costs.push_back(banded_trim(a,b));
    std::vector<Result> traces{hirschberg(a,b,true), sequential(a,b,true), nvpd_omp(a,b,true)};
#ifdef CHPS_KOKKOS
    costs.push_back(nvpd_kokkos(a,b));
    traces.push_back(nvpd_kokkos(a,b,true));
#endif
#ifdef CHPS_CUDA
    costs.push_back(nvpd_cuda(a,b));
    traces.push_back(nvpd_cuda(a,b,true));
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
