#pragma once

#include "corridor.cuh"
#include <cooperative_groups.h>
#include <cstdint>

namespace chps::detail::corridor {

// One cooperative block can own several column tiles when the full band
// exceeds the number of resident blocks. All blocks join the same grid barrier.
__global__ void persistent_tile(int* frontier, const Chunk* index, const int* edges,
                                const std::uint32_t* query, int m, int n, int margin, int right,
                                int upper, int capacity, int pitch, int first_slab, int slabs) {
    const auto grid = cooperative_groups::this_grid();
    extern __shared__ int scratch[];
    Chunk* cache = reinterpret_cast<Chunk*>(scratch + 2 * pitch);
    for (int slab = first_slab; slab < first_slab + slabs; ++slab) {
        const int first = slab * rows + 1;
        const int depth = min(rows, m - first + 1);
        const int end_row = first + depth - 1;
        const int lo = max(0, end_row - margin) / columns * columns;
        const int last_col_tile = min(n, end_row + right) / columns;
        const int prev_lo = max(0, first - 1 - margin) / columns * columns;
        const int prev_hi = slab == 0 ? min(n, right) :
            min(n, (min(n, first - 1 + right) / columns + 1) * columns - 1);
        const int* input = slab == 0 ? nullptr : frontier + ((slab - 1) & 1) * capacity;
        int* output = frontier + (slab & 1) * capacity;
        for (int tile_id = lo / columns + blockIdx.x; tile_id <= last_col_tile; tile_id += gridDim.x) {
            const int begin = tile_id * columns;
            const int end = min(n, begin + columns - 1);
            const int* left = edges + static_cast<std::size_t>(tile_id) * rows + (rows - depth);
            const int first_chunk = left[0] / 32;
            for (int c = threadIdx.x; c <= end / 32 - first_chunk; c += blockDim.x)
                cache[c] = index[first_chunk + c];
            __syncthreads();
            const int* previous = input;
            int previous_start = prev_lo;
            for (int t = 0; t < depth; ++t) {
                const int row = first + t;
                const int base = row - 1;
                const int symbol = (query[base / 16] >> (2 * (base & 15))) & 3u;
                const int origin = left[t];
                int* next = scratch + (t & 1) * pitch;
                for (int x = threadIdx.x; x <= end - origin; x += blockDim.x) {
                    const int j = origin + x;
                    int value = row;
                    if (j) {
                        const int k = last(cache[j / 32 - first_chunk], j, symbol);
                        if (t) {
                            value = nvpd_cell(previous, j - previous_start, k ? k - previous_start : 0);
                        } else {
                            const int diag = exterior(input, j-1, prev_lo, prev_hi, upper+1);
                            value = diag;
                            if (k != j) {
                                value = 1 + min(diag, exterior(input, j, prev_lo, prev_hi, upper+1));
                                if (k) value = min(value, exterior(input, k-1, prev_lo, prev_hi, upper+1) + (j-k));
                            }
                        }
                    }
                    if (t + 1 == depth) output[j - lo] = value;
                    else next[x] = value;
                }
                __syncthreads();
                previous = next;
                previous_start = origin;
            }
        }
        // At the batch boundary the next launch, on the same stream, supplies
        // this dependency. Inside a batch every CTA must publish its frontier.
        if (slab + 1 < first_slab + slabs) grid.sync();
    }
}

template<class Buffer>
std::optional<Result> persistent_score(std::string_view a, std::string_view b, Buffer& index,
                                       Buffer& table, Buffer& boundaries) {
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const auto layout = make_layout(a, b);
    if (layout.path == Path::length_difference) return Result{layout.upper, {}, 0};
    if (layout.path != Path::tiles) return std::nullopt;
    int device, cooperative;
    check(cudaGetDevice(&device), "get persistent CUDA device");
    check(cudaDeviceGetAttribute(&cooperative, cudaDevAttrCooperativeLaunch, device),
          "check cooperative launch support");
    if (!cooperative) return std::nullopt;
    const int pitch = layout.pitch, capacity = layout.capacity;
    const std::size_t smem = 2*pitch*sizeof(int) + ((pitch+62)/32)*sizeof(Chunk);
    int blocks_per_sm, sms;
    check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&blocks_per_sm, persistent_tile, threads, smem),
          "get cooperative occupancy");
    check(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device), "count CUDA SMs");
    const int blocks = std::min((capacity + columns - 1) / columns, sms * blocks_per_sm);
    if (blocks == 0) return std::nullopt;

    std::vector<std::uint32_t> packed((a.size()+15)/16, 0);
    for (int i = 0; i < m; ++i)
        packed[i/16] |= static_cast<std::uint32_t>(dna_id(a[i])) << (2*(i&15));
    index.resize(layout.index.size()*sizeof(Chunk)/sizeof(int));
    table.resize(2*static_cast<std::size_t>(capacity));
    boundaries.resize(layout.edges.size() + packed.size());
    auto* device_query = reinterpret_cast<std::uint32_t*>(boundaries.data + layout.edges.size());
    check(cudaMemcpy(index.data, layout.index.data(), index.bytes, cudaMemcpyHostToDevice), "copy persistent index");
    check(cudaMemcpy(boundaries.data, layout.edges.data(), layout.edges.size()*sizeof(int), cudaMemcpyHostToDevice),
          "copy persistent halos");
    check(cudaMemcpy(device_query, packed.data(), packed.size()*sizeof(std::uint32_t), cudaMemcpyHostToDevice),
          "copy persistent query");
    // Cooperative launches take pointers to storage holding each argument.
    // The runtime copies those arguments synchronously during the launch.
    const Chunk* device_index = reinterpret_cast<const Chunk*>(index.data);
    int upper = layout.upper, margin = (upper-(n-m))/2, right = n-m+margin;
    int length_a = m, length_b = n, row_capacity = capacity, shared_pitch = pitch;
    const int total_slabs = (m+rows-1)/rows;
    constexpr int batch_slabs = 256;
    for (int slab = 0; slab < total_slabs;) {
        int count = std::min(batch_slabs,total_slabs-slab);
        void* args[] = {&table.data, &device_index, &boundaries.data, &device_query, &length_a, &length_b,
                        &margin, &right, &upper, &row_capacity, &shared_pitch, &slab, &count};
        const auto status = cudaLaunchCooperativeKernel(reinterpret_cast<const void*>(persistent_tile),
                                                         blocks, threads, args, smem, nullptr);
        if (slab == 0 && (status == cudaErrorCooperativeLaunchTooLarge || status == cudaErrorNotSupported)) {
            cudaGetLastError();
            return std::nullopt;
        }
        check(status, "launch persistent NvPD corridor");
        slab += count;
    }
    const auto slice = output_slice(m,m,n,upper);
    const int* terminal = table.data + ((total_slabs-1)&1)*capacity + (n-slice.first);
    int score;
    check(cudaMemcpy(&score,terminal,sizeof(score),cudaMemcpyDeviceToHost), "read persistent score");
    if (score > upper) throw std::logic_error("persistent score exceeded its constructive bound");
    return Result{score,{},index.bytes+table.bytes+boundaries.bytes};
}
} // namespace chps::detail::corridor
