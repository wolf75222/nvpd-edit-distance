#pragma once

#include "local_tiles.cuh"
#include "packed_index.hpp"

namespace chps::detail {
template<class Symbols>
__global__ void packed_local_tile(const int* input, int* output, const int* checkpoints,
                                  const int* starts, const unsigned* words, Symbols symbols,
                                  int first_row, int count, int n, int columns, int pitch) {
    extern __shared__ int scores[];
    unsigned* cached = reinterpret_cast<unsigned*>(scores + 2 * pitch);
    const int begin = static_cast<int>(blockIdx.x) * columns;
    const int end = begin + min(n - begin, columns - 1);
    const int* edges = starts + blockIdx.x * local_tile_rows + local_tile_rows - count;
    const int origin = edges[0];
    for (int x = threadIdx.x; x <= end - origin; x += blockDim.x)
        cached[x] = words[origin + x];
    __syncthreads();
    const int* prev = input;
    int previous_start = 0;
    for (int t = 0; t < count; ++t) {
        const int left = edges[t], symbol = symbols.id[t];
        int* curr = scores + (t % 2) * pitch;
        for (int x = threadIdx.x; x <= end - left; x += blockDim.x) {
            const int j = left + x;
            const int k = unpack_match(cached[j - origin], checkpoints, j, symbol);
            const int value = j ? nvpd_cell(prev, j - previous_start, k ? k - previous_start : 0)
                                : first_row + t;
            if (t + 1 == count) output[j] = value;
            else curr[x] = value;
        }
        __syncthreads();
        prev = curr;
        previous_start = left;
    }
}
} // namespace chps::detail
