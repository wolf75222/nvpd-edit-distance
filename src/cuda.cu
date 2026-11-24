#include <cuda_runtime.h>

#include "dp_internal.hpp"

#include <algorithm>
#include <stdexcept>

namespace chps {
namespace {

void check(cudaError_t status, const char* action) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(action) + ": " + cudaGetErrorString(status));
}

struct DeviceInts {
    int* data = nullptr;
    const std::size_t bytes;
    explicit DeviceInts(std::size_t count) : bytes(count * sizeof(int)) {
        check(cudaMalloc(&data, bytes), "cudaMalloc");
    }
    ~DeviceInts() { cudaFree(data); }
    DeviceInts(const DeviceInts&) = delete;
    DeviceInts& operator=(const DeviceInts&) = delete;
};

__global__ void first_row(int* row, int n) {
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j <= n) row[j] = j;
}

__global__ void next_row(const int* prev, int* curr, const int* last, int i, int n) {
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j <= n) curr[j] = j == 0 ? i : detail::nvpd_cell(prev, j, last[j]);
}

} // namespace

Result nvpd_cuda(std::string_view a, std::string_view b, bool trace) {
    detail::validate(a, b, trace);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    if (!trace && a.size() > b.size()) std::swap(a, b);
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const std::size_t stride = b.size() + 1;
    const detail::MatchIndex matches(a, b);
    DeviceInts indices(matches.entries.size());
    check(cudaMemcpy(indices.data, matches.entries.data(), indices.bytes,
                     cudaMemcpyHostToDevice), "copy match index");
    const std::size_t table_size = (trace ? a.size() + 1 : 2) * stride;
    DeviceInts table(table_size);
    constexpr int threads = 256;
    const int blocks = static_cast<int>((stride + threads - 1) / threads);
    first_row<<<blocks, threads>>>(table.data, n);
    check(cudaGetLastError(), "launch first row");
    int* prev = table.data;
    // Launches use one stream, so each row sees the completed preceding row.
    for (int i = 1; i <= m; ++i) {
        int* curr = table.data + (trace ? i : i % 2) * stride;
        const auto s = matches.symbol[static_cast<unsigned char>(a[i-1])];
        next_row<<<blocks, threads>>>(prev, curr, indices.data + s * stride, i, n);
        check(cudaGetLastError(), "launch NvPD row");
        prev = curr;
    }
    const auto device_bytes = indices.bytes + table.bytes;
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

} // namespace chps
