#pragma once

#include "dp_internal.hpp"
#include "corridor_plan.hpp"
#include <cuda_runtime.h>
#include <optional>
#include <stdexcept>
#include <string>

namespace chps::detail {
namespace corridor {

struct Slab {
    unsigned char symbol[rows];
};

__device__ inline int exterior(const int* input, int j, int lo, int hi, int inf) {
    if (j < lo || j > hi) return inf;
    return input ? input[j-lo] : j;
}

// The first row receives a possibly truncated previous frontier. Each later
// row has all its predecessors in the halo stored in shared memory.
__global__ void tile(const int* input, int* output, const Chunk* index, const int* edges,
                     Slab slab, int first, int count, int n, int in_lo, int in_hi,
                     int out_lo, int pitch, int inf) {
    extern __shared__ int smem[];
    const int begin = out_lo + static_cast<int>(blockIdx.x) * columns;
    const int end = min(n,begin+columns-1);
    const int* edge = edges + static_cast<std::size_t>(begin / columns) * rows + (rows-count);
    const int first_chunk = edge[0] / 32;
    Chunk* cache = reinterpret_cast<Chunk*>(smem + 2*pitch);
    for (int c = threadIdx.x; c <= end/32-first_chunk; c += blockDim.x)
        cache[c] = index[first_chunk+c];
    __syncthreads();
    const int* previous = input;
    int prev_lo = in_lo;
    for (int t = 0; t < count; ++t) {
        const int lo = edge[t];
        int* next = smem + (t&1)*pitch;
        for (int x = threadIdx.x; x <= end-lo; x += blockDim.x) {
            const int j = lo+x;
            int v = first+t;
            if (j) {
                const int k = last(cache[j/32-first_chunk],j,slab.symbol[t]);
                if (t) {
                    v = nvpd_cell(previous,j-prev_lo,k ? k-prev_lo : 0);
                } else {
                    const int diag = exterior(previous,j-1,in_lo,in_hi,inf);
                    if (k == j) v = diag;
                    else {
                        v = 1 + min(diag,exterior(previous,j,in_lo,in_hi,inf));
                        if (k) v = min(v,exterior(previous,k-1,in_lo,in_hi,inf) + (j-k));
                    }
                }
            }
            if (t+1 == count) output[j-out_lo] = v;
            else next[x] = v;
        }
        __syncthreads();
        previous = next;
        prev_lo = lo;
    }
}

inline void check(cudaError_t err, const char* action) {
    if (err != cudaSuccess)
        throw std::runtime_error(std::string(action) + ": " + cudaGetErrorString(err));
}
} // namespace corridor

// Any path of cost U crossing (i,j) pays at least |j-i| insertions/deletions
// before the crossing and |(n-m)-(j-i)| afterwards. Rounding each output
// interval outwards to full tiles therefore retains at least one optimal path.
template<class Buffer>
std::optional<Result> nvpd_corridor(std::string_view a, std::string_view b, Buffer& index,
                                    Buffer& table, Buffer& boundaries, std::size_t peak) {
    using namespace corridor;
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const Layout layout = make_layout(a,b);
    if (layout.path == Path::length_difference) return Result{layout.upper,{},peak};
    if (layout.path != Path::tiles) return std::nullopt;
    const int capacity = layout.capacity, upper = layout.upper, pitch = layout.pitch;
    const auto& chunks = layout.index;
    const auto& starts = layout.edges;
    auto resident = [&]{return index.bytes + table.bytes + boundaries.bytes;};
    index.resize(chunks.size()*sizeof(Chunk)/sizeof(int));
    peak = std::max(peak,resident());
    table.resize(2*static_cast<std::size_t>(capacity));
    peak = std::max(peak,resident());
    boundaries.resize(starts.size());
    peak = std::max(peak,resident());
    corridor::check(cudaMemcpy(index.data,chunks.data(),index.bytes,cudaMemcpyHostToDevice),"copy corridor index");
    corridor::check(cudaMemcpy(boundaries.data,starts.data(),boundaries.bytes,cudaMemcpyHostToDevice),"copy corridor halos");
    int in_lo = 0, in_hi = n-m+(upper-(n-m))/2;
    int* prev = nullptr;
    int* next = table.data;
    const std::size_t smem = 2*pitch*sizeof(int) + ((pitch+62)/32)*sizeof(Chunk);
    for (int i = 1; i <= m; ) {
        const int count = std::min(rows,m-i+1);
        const int end_row = i+count-1;
        const auto [lo,hi] = output_slice(end_row,m,n,upper);
        Slab slab{};
        for (int t = 0; t < count; ++t) slab.symbol[t] = static_cast<unsigned char>(dna_id(a[i+t-1]));
        tile<<<hi/columns-lo/columns+1,threads,smem>>>(prev,next,reinterpret_cast<const Chunk*>(index.data),
            boundaries.data,slab,i,count,n,in_lo,in_hi,lo,pitch,upper+1);
        corridor::check(cudaGetLastError(),"launch NvPD corridor");
        prev = next;
        next = prev == table.data ? table.data+capacity : table.data;
        in_lo = lo;
        in_hi = hi;
        i += count;
    }
    int score;
    corridor::check(cudaMemcpy(&score,prev+(n-in_lo),sizeof(score),cudaMemcpyDeviceToHost),"copy corridor score");
    if (score > upper) throw std::logic_error("NvPD corridor exceeded its constructive bound");
    return Result{score,{},peak};
}
} // namespace chps::detail
