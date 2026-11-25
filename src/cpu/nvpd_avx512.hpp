#pragma once

#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__)) && !defined(__CUDACC__)
#define CHPS_NVPD_AVX512_DISPATCH 1

#include "chps/dp.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <omp.h>
#include <stdexcept>
#include <string_view>
#include <vector>
#include <immintrin.h>

namespace chps::detail::avx512 {
inline constexpr int block_lanes = 16;

struct MatchBlock {
    std::array<std::uint8_t, block_lanes> source_lane;
    int preceding_match;
};

__attribute__((target("avx512f"))) inline void evaluate_block(int* current, const int* diagonal, int carry,
                           const MatchBlock& block) {
    const auto one = _mm512_set1_epi32(1);
    const auto diag = _mm512_loadu_si512(diagonal);
    const auto top = _mm512_loadu_si512(diagonal + 1);
    const auto labels = _mm512_maskz_cvtepu8_epi32(0xffff, _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(block.source_lane.data())));
    const auto last = _mm512_maskz_permutex2var_epi32(0xffff, diag, labels, _mm512_set1_epi32(carry));
    const auto ordinary = _mm512_maskz_min_epi32(0xffff, diag, _mm512_add_epi32(top, one));
    _mm512_storeu_si512(current, _mm512_maskz_min_epi32(0xffff, ordinary, _mm512_sub_epi32(last, one)));
}

__attribute__((target("avx512f"))) inline void evaluate_span(int* curr, const int* prev,
        const MatchBlock* matches, std::size_t begin, std::size_t end) {
    for (auto chunk = begin; chunk < end; ++chunk) {
        const auto first = chunk * block_lanes + 1;
        const auto& block = matches[chunk];
        evaluate_block(curr + first + 1, prev + first, prev[block.preceding_match], block);
    }
}

// Row position j+1 holds F(i,j) = D(i,j) - j. The last match contributes
// F(i-1,k-1) - 1. If k=0, slot zero contributes i-1, an upper bound for j>=1.
// Inside each block its source is either a local lane or the single match
// preceding the whole block. The padded final block writes into the same row.
inline Result score(std::string_view a, std::string_view b) {
    const auto blocks = (b.size() + block_lanes - 1) / block_lanes;
    std::array<unsigned char, 256> alphabet;
    std::array<int, 256> symbols;
    symbols.fill(-1);
    int count = 0;
    for (unsigned char c : a) {
        if (symbols[c] < 0) {
            alphabet[count] = c;
            symbols[c] = count++;
        }
    }
    std::vector<MatchBlock> index;
    if (count && blocks > index.max_size() / static_cast<std::size_t>(count))
        throw std::length_error("blocked match index exceeds addressable memory");
    index.resize(blocks * static_cast<std::size_t>(count));
    #pragma omp parallel for schedule(static)
    for (int s = 0; s < count; ++s) {
        int last = 0;
        for (std::size_t chunk = 0; chunk < blocks; ++chunk) {
            const auto offset = chunk * block_lanes;
            auto& block = index[s * blocks + chunk];
            block.preceding_match = last;
            int source = block_lanes;  // The second permutation operand broadcasts the preceding match.
            for (int lane = 0; lane < block_lanes; ++lane) {
                const auto position = offset + lane;
                if (position < b.size() && static_cast<unsigned char>(b[position]) == alphabet[s]) {
                    last = static_cast<int>(position + 1);
                    source = lane;
                }
                block.source_lane[lane] = static_cast<std::uint8_t>(source);
            }
        }
    }
    const std::size_t pitch = blocks * block_lanes + 2;
    std::vector<int> rows;
    if (pitch > rows.max_size() / 2) throw std::length_error("score rows exceed addressable memory");
    rows.resize(2 * pitch);
    rows[0] = 1;
    #pragma omp parallel
    {
        const auto workers = static_cast<std::size_t>(omp_get_num_threads());
        const auto worker = static_cast<std::size_t>(omp_get_thread_num());
        const auto quotient = blocks / workers, extra = blocks % workers;
        const auto begin = worker * quotient + std::min(worker, extra);
        const auto end = begin + quotient + (worker < extra);
        for (int i = 1; i <= static_cast<int>(a.size()); ++i) {
            const int* prev = rows.data() + ((i - 1) % 2) * pitch;
            int* curr = rows.data() + (i % 2) * pitch;
            if (worker == 0) {
                curr[0] = i + 1;
                curr[1] = i;
            }
            const auto s = symbols[static_cast<unsigned char>(a[i-1])];
            const auto* matches = index.data() + s * blocks;
            evaluate_span(curr, prev, matches, begin, end);
            #pragma omp barrier
        }
    }
    return {rows[(a.size() % 2) * pitch + b.size() + 1] + static_cast<int>(b.size()), "", 0};
}
} // namespace chps::detail::avx512
#endif
