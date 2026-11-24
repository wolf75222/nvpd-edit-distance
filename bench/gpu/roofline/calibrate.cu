#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cuda_runtime.h>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using U = std::uint32_t;
constexpr int threads = 256;
constexpr int words = 1 << 25;

static void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}

// Interior NvPD work unit: two minima, two additions, one subtraction.
// Volatile PTX retains every update; SASS may combine the integer additions.
__device__ __forceinline__ U update(U up, U diagonal, U match, U column, U previous) {
    U result;
    asm volatile("{ .reg .u32 a,b; min.u32 a,%1,%2; add.u32 a,a,1; "
                 "sub.u32 b,%4,%5; add.u32 b,b,%3; min.u32 %0,a,b; }"
                 : "=r"(result)
                 : "r"(up), "r"(diagonal), "r"(match), "r"(column), "r"(previous));
    return result;
}

__global__ void fill(U *a, U *b, U *c) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    a[i] = i * 17u + 1;
    b[i] = i * 7u + 3;
    c[i] = i * 5u + 11;
}

__global__ void stream_copy(const uint4 *a, uint4 *out) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    out[i] = a[i];
}

template <int N> __global__ void stream_updates(const U *a, const U *b, const U *c, U *out) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    U x = a[i], y = b[i], z = c[i], v = x;
#pragma unroll 4
    for (int k = 0; k < N; ++k)
        v = update(v, y, z, v, x);
    out[i] = v;
}

__global__ void register_updates(U *out) {
    const U i = blockIdx.x * blockDim.x + threadIdx.x;
    U v[8];
#pragma unroll
    for (int j = 0; j < 8; ++j)
        v[j] = i * (j + 5) + 1;
#pragma unroll 4
    for (int r = 0; r < 512; ++r) {
#pragma unroll
        for (int j = 0; j < 8; ++j)
            v[j] = update(v[j], i + j, i * 7 + j, v[j], i * 17 + j);
    }
    U sum = 0;
#pragma unroll
    for (int j = 0; j < 8; ++j)
        sum += v[j];
    out[i] = sum;
}

template <class Launch>
void time_kernel(const char *name, int block, int steps, unsigned long long updates,
                 unsigned long long bytes, Launch launch) {
    cudaEvent_t start, stop;
    check(cudaEventCreate(&start));
    check(cudaEventCreate(&stop));
    launch();
    check(cudaDeviceSynchronize());
    constexpr int repeats = 7;
    check(cudaEventRecord(start));
    for (int k = 0; k < repeats; ++k)
        launch();
    check(cudaEventRecord(stop));
    check(cudaEventSynchronize(stop));
    check(cudaGetLastError());
    float ms;
    check(cudaEventElapsedTime(&ms, start, stop));
    std::cout << name << ',' << block << ',' << steps << ',' << updates << ',' << bytes << ','
              << std::setprecision(12) << ms / (1000.0 * repeats) << '\n';
    check(cudaEventDestroy(start));
    check(cudaEventDestroy(stop));
}

template <int N> void stream(int block, U *a, U *b, U *c, U *dst) {
    time_kernel("stream-updates", block, N, 1ull * words * N, 16ull * words,
                [&] { stream_updates<N><<<words / threads, threads>>>(a, b, c, dst); });
    // Independently check a first, interior and last word after every sweep.
    for (U i : {0u, 12345u, static_cast<U>(words - 1)}) {
        U x = i * 17 + 1, y = i * 7 + 3, z = i * 5 + 11, expected = x, observed;
        for (int k = 0; k < N; ++k)
            expected = std::min(std::min(expected, y) + 1u, z + (expected - x));
        check(cudaMemcpy(&observed, dst + i, sizeof(U), cudaMemcpyDeviceToHost));
        if (observed != expected)
            throw std::runtime_error("stream checksum mismatch");
    }
}

int main(int argc, char **argv) try {
    if (argc != 2)
        return 2;
    const int block = std::stoi(argv[1]);
    U *a, *b, *c, *dst;
    for (U **p : {&a, &b, &c, &dst})
        check(cudaMalloc(p, words * sizeof(U)));
    fill<<<words / threads, threads>>>(a, b, c);
    check(cudaDeviceSynchronize());
    time_kernel("stream-copy", block, 0, 0, 8ull * words, [&] {
        stream_copy<<<words / (4 * threads), threads>>>(reinterpret_cast<uint4 *>(a),
                                                        reinterpret_cast<uint4 *>(dst));
    });
    constexpr int reg_threads = 1 << 18;
    time_kernel("register-updates", block, 4096, 4096ull * reg_threads, 4ull * reg_threads,
                [&] { register_updates<<<reg_threads / threads, threads>>>(dst); });
    stream<1>(block, a, b, c, dst);
    stream<4>(block, a, b, c, dst);
    stream<16>(block, a, b, c, dst);
    stream<64>(block, a, b, c, dst);
    stream<256>(block, a, b, c, dst);
    for (U *p : {a, b, c, dst})
        check(cudaFree(p));
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
