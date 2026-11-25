#pragma once

#include "dp_internal.hpp"
#include <algorithm>

namespace chps::detail {
inline std::vector<unsigned> pack_match_offsets(const MatchIndex& index, int n) {
    const std::size_t stride = static_cast<std::size_t>(n) + 1;
    std::vector<unsigned> words(stride);
    for (unsigned sym = 0; sym < 4; ++sym) {
        const int* last = index.entries.data() + sym * stride;
        for (int j = 0; j <= n; ++j) {
            const unsigned delta = last[j] ? std::min(j - last[j], 255) : 255;
            words[j] |= delta << (8 * sym);
        }
    }
    return words;
}

inline constexpr int checkpoint_stride = 256;

// For j in [256b, 256b+255], j-k >= 255 implies k <= 256b. The
// last occurrence at 256b resolves all escapes, including absent symbols.
inline std::vector<int> match_checkpoints(const MatchIndex& index, int n) {
    const std::size_t stride = static_cast<std::size_t>(n) + 1;
    const std::size_t blocks = (stride + checkpoint_stride - 1) / checkpoint_stride;
    std::vector<int> checkpoints(4 * blocks);
    for (std::size_t block = 0; block < blocks; ++block)
        for (std::size_t symbol = 0; symbol < 4; ++symbol)
            checkpoints[4 * block + symbol] = index.entries[symbol * stride + block * checkpoint_stride];
    return checkpoints;
}

#ifdef __CUDACC__
__host__ __device__
#endif
inline int unpack_match(unsigned word, const int* checkpoints, int j, int symbol) {
    const unsigned delta = (word >> (8 * symbol)) & 255u;
    return delta == 255u ? checkpoints[4 * (j / checkpoint_stride) + symbol]
                         : j - static_cast<int>(delta);
}
} // namespace chps::detail
