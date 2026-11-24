#pragma once

#include "dp_internal.hpp"
#include <cuda_runtime.h>
#include <algorithm>

namespace chps::detail {
#ifndef CHPS_CUDA_TILE_ROWS
#define CHPS_CUDA_TILE_ROWS 8
#endif
inline constexpr int local_tile_rows = CHPS_CUDA_TILE_ROWS;
static_assert(local_tile_rows >= 1 && local_tile_rows <= 32);

// Row r is produced for [left[r], end] in each block. k_s(j), the last
// occurrence of symbol s at or before j, is nondecreasing once it exists.
// Thus at the left edge L, min(L-1, min_s(k_s(L)-1 for k_s(L)>0)) bounds all
// predecessor columns for j >= L. If k_s(L)==0, any later match lies right of L.
struct LocalTilePlan {
    std::vector<int> left;
    int pitch = 0;
};

inline LocalTilePlan make_local_plan(const MatchIndex& index, int n, int columns, int maximum_pitch) {
    const std::size_t stride = static_cast<std::size_t>(n) + 1;
    const std::size_t blocks = (stride + columns - 1) / columns;
    LocalTilePlan plan{std::vector<int>(blocks * local_tile_rows)};
    for (std::size_t block = 0; block < blocks; ++block) {
        int* left = plan.left.data() + block * local_tile_rows;
        const int begin = static_cast<int>(block * columns);
        left[local_tile_rows-1] = begin;
        for (int t = local_tile_rows - 1; t > 0; --t) {
            const int column = left[t];
            int predecessor = std::max(0, column - 1);
            for (std::size_t row = 0; row < index.entries.size(); row += stride) {
                const int k = index.entries[row + column];
                if (k) predecessor = std::min(predecessor, k - 1);
            }
            left[t-1] = predecessor;
        }
        const int out_width = std::min(columns, n - begin + 1);
        plan.pitch = std::max(plan.pitch, begin - left[0] + out_width);
        if (plan.pitch > maximum_pitch) return {};
    }
    return plan;
}

template<class Symbols>
__global__ void local_row_tile(const int* input, int* output, const int* indices, const int* starts,
                               Symbols symbols, int first_row, int count, int n, int columns, int pitch) {
    extern __shared__ int scores[];
    const int begin = static_cast<int>(blockIdx.x) * columns;
    const int end = begin + min(n - begin, columns - 1);
    const int* left_edges = starts + blockIdx.x * local_tile_rows + local_tile_rows - count;
    const int* prev = input;
    int prev_start = 0;
    for (int t = 0; t < count; ++t) {
        const int left = left_edges[t];
        const int* last = indices + static_cast<std::size_t>(symbols.id[t]) * (n + 1);
        int* curr = scores + (t % 2) * pitch;
        for (int x = threadIdx.x; x <= end - left; x += blockDim.x) {
            const int j = left + x;
            const int k = last[j];
            const int value = j ? nvpd_cell(prev, j - prev_start, k ? k - prev_start : 0) : first_row + t;
            if (t + 1 == count) output[j] = value;
            else curr[j - left] = value;
        }
        __syncthreads();
        prev = curr;
        prev_start = left;
    }
}
} // namespace chps::detail
