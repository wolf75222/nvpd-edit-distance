#include <cuda_runtime.h>

#include "chps/cuda.hpp"
#include "dp_internal.hpp"
#include "local_tiles.cuh"
#include "packed_match.cuh"

#ifndef CHPS_CUDA_PACKED_MATCHES
#define CHPS_CUDA_PACKED_MATCHES 0
#endif

#include <algorithm>
#include <stdexcept>

namespace chps {
namespace {

void check(cudaError_t status, const char* action) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(action) + ": " + cudaGetErrorString(status));
}

__global__ void first_row(int* row, int n) {
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j <= n) row[j] = j;
}

__global__ void next_row(const int* prev, int* curr, const int* last, int i, int n) {
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j <= n) curr[j] = j == 0 ? i : detail::nvpd_cell(prev, j, last[j]);
}

constexpr int tile_rows = detail::local_tile_rows;
constexpr int max_halo_step = 64;
struct Symbols { int id[tile_rows]; };

// Largest leftward dependency in any MIT row. A missing match uses j and j-1.
int left_reach(const detail::MatchIndex& matches, std::size_t stride) {
    int reach = 1;
    for (std::size_t offset = 0; offset < matches.entries.size(); offset += stride) {
        for (int j = 1; j < static_cast<int>(stride); ++j) {
            const int k = matches.entries[offset + j];
            if (k) reach = std::max(reach, j - k + 1);
            if (reach > max_halo_step) return reach;
        }
    }
    return reach;
}

// After t rows, column j can only depend on columns [j-t*reach, j].
// Redundant cells at the left of each block let it advance without a grid barrier.
__global__ void row_tile(const int* input, int* output, const int* indices,
                         Symbols symbols, int first_i, int count, int n, int reach,
                         int output_width) {
    extern __shared__ int scratch[];
    const int width = static_cast<int>(blockDim.x);
    const int begin = static_cast<int>(blockIdx.x) * output_width;
    const int end = begin + min(n - begin, output_width - 1);
    const int pitch = output_width + (count - 1) * reach;
    const int* previous = input;
    int previous_begin = 0;
    for (int t = 0; t < count; ++t) {
        const int left = max(0, begin - (count - t - 1) * reach);
        const int* last = indices + static_cast<std::size_t>(symbols.id[t]) * (n + 1);
        int* next = scratch + (t % 2) * pitch;
        for (int x = threadIdx.x; x <= end - left; x += width) {
            const int j = x + left;
            const int k = last[j];
            const int value = j == 0 ? first_i + t :
                detail::nvpd_cell(previous, j - previous_begin, k ? k - previous_begin : 0);
            if (t + 1 == count) output[j] = value;
            else next[x] = value;
        }
        __syncthreads();
        previous = next;
        previous_begin = left;
    }
}

} // namespace

void CudaWorkspace::Buffer::resize(std::size_t count) {
    const auto requested_bytes = count * sizeof(int);
    if (bytes == requested_bytes) return;
    if (data) {
        check(cudaFree(data), "release previous CUDA buffer");
        data = nullptr;
        bytes = 0;
    }
    int* replacement = nullptr;
    check(cudaMalloc(&replacement, requested_bytes), "allocate CUDA buffer");
    data = replacement;
    bytes = requested_bytes;
}

CudaWorkspace::Buffer::~Buffer() {
    if (data) cudaFree(data);
}

Result CudaWorkspace::distance(std::string_view a, std::string_view b, bool trace) {
    detail::validate(a, b, trace);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    if (!trace && a.size() > b.size()) std::swap(a, b);
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const std::size_t stride = b.size() + 1;
    const detail::MatchIndex matches(a, b);
    auto& indices = indices_;
    indices.resize(matches.entries.size());
    check(cudaMemcpy(indices.data, matches.entries.data(), indices.bytes,
                     cudaMemcpyHostToDevice), "copy match index");
    const std::size_t table_size = (trace ? a.size() + 1 : 2) * stride;
    auto& table = table_;
    table.resize(table_size);
    constexpr int threads = 256;
    const int blocks = static_cast<int>((stride + threads - 1) / threads);
    first_row<<<blocks, threads>>>(table.data, n);
    check(cudaGetLastError(), "launch first row");
    int* prev = table.data;
    const int reach = !trace && m >= tile_rows ? left_reach(matches, stride) : max_halo_step + 1;
    const bool consider_local = halo_ == CudaHalo::local && !trace && m >= tile_rows && reach > max_halo_step;
    int columns = threads;
    detail::LocalTilePlan plan;
    bool pack = false;
    if (reach <= max_halo_step || consider_local) {
        int device, multiprocessors;
        check(cudaGetDevice(&device), "get CUDA device");
        check(cudaDeviceGetAttribute(&multiprocessors, cudaDevAttrMultiProcessorCount, device),
              "get SM count");
        // Wider blocks amortize the halo when at least two CTAs per SM remain.
        const int wide = 4 * threads;
        columns = stride / wide >= 2 * multiprocessors ? wide : threads;
        if (consider_local) {
            int shared_limit;
            check(cudaDeviceGetAttribute(&shared_limit, cudaDevAttrMaxSharedMemoryPerBlock, device),
                  "get CUDA shared memory limit");
            plan = detail::make_local_plan(matches, n, columns, shared_limit / (2 * sizeof(int)));
            pack = CHPS_CUDA_PACKED_MATCHES && !plan.left.empty() && matches.entries.size() == 4 * stride
                && plan.pitch <= shared_limit / (3 * static_cast<int>(sizeof(int)));
        }
    }
    if (!plan.left.empty()) {
        boundaries_.resize(plan.left.size() + (pack ? stride : 0));
        if (pack) {
            auto words = detail::pack_match_offsets(matches, n);
            check(cudaMemcpy(boundaries_.data + plan.left.size(), words.data(), stride * sizeof(unsigned),
                             cudaMemcpyHostToDevice), "copy compact match offsets");
        }
        check(cudaMemcpy(boundaries_.data, plan.left.data(), plan.left.size() * sizeof(int), cudaMemcpyHostToDevice),
              "copy local tile boundaries");
    }
    if (reach <= max_halo_step || !plan.left.empty()) {
        const int tile_blocks = static_cast<int>((stride + columns - 1) / columns);
        int* curr = table.data + stride;
        for (int i = 1; i <= m;) {
            const int count = std::min(tile_rows, m - i + 1);
            Symbols symbols{};
            for (int t = 0; t < count; ++t)
                symbols.id[t] = matches.symbol[static_cast<unsigned char>(a[i-1+t])];
            if (plan.left.empty()) {
                const auto bytes = 2 * (columns + (count - 1) * reach) * sizeof(int);
                row_tile<<<tile_blocks, threads, bytes>>>(prev, curr, indices.data, symbols,
                                                          i, count, n, reach, columns);
            } else if (pack) {
                const auto bytes = 3 * plan.pitch * sizeof(int);
                const auto words = reinterpret_cast<const unsigned*>(boundaries_.data + plan.left.size());
                detail::packed_local_tile<<<tile_blocks, threads, bytes>>>(prev, curr, indices.data,
                    boundaries_.data, words, symbols, i, count, n, columns, plan.pitch);
            } else {
                const auto bytes = 2 * plan.pitch * sizeof(int);
                detail::local_row_tile<<<tile_blocks, threads, bytes>>>(prev, curr, indices.data,
                    boundaries_.data, symbols, i, count, n, columns, plan.pitch);
            }
            check(cudaGetLastError(), "launch NvPD tile");
            std::swap(prev, curr);
            i += count;
        }
    } else {
        for (int i = 1; i <= m; ++i) {
            int* curr = table.data + (trace ? i : i % 2) * stride;
            const auto s = matches.symbol[static_cast<unsigned char>(a[i-1])];
            next_row<<<blocks, threads>>>(prev, curr, indices.data + s * stride, i, n);
            check(cudaGetLastError(), "launch NvPD row");
            prev = curr;
        }
    }
    const auto device_bytes = indices.bytes + table.bytes + boundaries_.bytes;
    if (trace) {
        std::vector<int> host(table_size);
        check(cudaMemcpy(host.data(), table.data, table.bytes, cudaMemcpyDeviceToHost),
              "copy traceback matrix");
        return {host.back(), detail::traceback(host.data(), a, b), device_bytes};
    }
    int cost;
    check(cudaMemcpy(&cost, prev + n, sizeof(cost), cudaMemcpyDeviceToHost), "copy distance");
    return {cost, "", device_bytes};
}

Result nvpd_cuda(std::string_view a, std::string_view b, bool trace) {
    CudaWorkspace workspace;
    return workspace.distance(a, b, trace);
}

Result nvpd_cuda_local(std::string_view a, std::string_view b, bool trace) {
    CudaWorkspace workspace(CudaHalo::local);
    return workspace.distance(a, b, trace);
}

} // namespace chps
