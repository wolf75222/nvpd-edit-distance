#pragma once

#include "dp_internal.hpp"
#include <cuda_runtime.h>

namespace chps::detail {

__device__ __forceinline__ int intermediate(const int *scores, const int *last, int row, int col) {
    return col ? nvpd_cell(scores, col, last[col]) : row;
}

// Recover the required entries of row i in registers and write only row i+1.
// Each dependency can reach any column of the previous globally completed row.
// Shuffles share intermediate values among the 32 consecutive lanes of a warp.
template <bool Shuffle>
__global__ void pair_rows(const int *prev, int *out, const int *first, const int *second, int i, int n) {
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if constexpr (!Shuffle) {
        if (j > n)
            return;
        if (!j) {
            out[0] = i + 1;
            return;
        }
        const int k = second[j];
        const int diagonal = intermediate(prev, first, i, j - 1);
        if (k == j) {
            out[j] = diagonal;
            return;
        }
        const int upper = intermediate(prev, first, i, j);
        const int edit = 1 + min(upper, diagonal);
        if (!k) {
            out[j] = edit;
            return;
        }
        out[j] = min(edit, intermediate(prev, first, i, k - 1) + (j - k));
    } else {
        // Invalid tail lanes remain active for the warp collectives and never
        // write. Any value requested by a valid lane comes from a valid lane.
        const int column = min(j, n);
        const int here = intermediate(prev, first, i, column);
        const int k = second[column];
        const int lane = threadIdx.x & 31;
        const int preceding = __shfl_up_sync(0xffffffff, here, 1);
        const int matched = __shfl_sync(0xffffffff, here, (k - 1) & 31);
        if (j > n)
            return;
        if (!j) {
            out[0] = i + 1;
            return;
        }
        const int diagonal = lane ? preceding : intermediate(prev, first, i, j - 1);
        if (k == j) {
            out[j] = diagonal;
            return;
        }
        const int edit = 1 + min(here, diagonal);
        if (!k) {
            out[j] = edit;
            return;
        }
        const int from_match = k - 1 >= j - lane ? matched : intermediate(prev, first, i, k - 1);
        out[j] = min(edit, from_match + (j - k));
    }
}

} // namespace chps::detail
