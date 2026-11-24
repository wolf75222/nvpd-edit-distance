#pragma once

#include "chps/dp.hpp"

namespace chps {

enum class CudaHalo { global, local };
Result nvpd_cuda_local(std::string_view a, std::string_view b, bool trace = false);

// Reuse exact-sized buffers over serial alignments on the same CUDA device.
class CudaWorkspace {
public:
    explicit CudaWorkspace(CudaHalo halo = CudaHalo::global) : halo_(halo) {}
    CudaWorkspace(const CudaWorkspace&) = delete;
    CudaWorkspace& operator=(const CudaWorkspace&) = delete;
    Result distance(std::string_view a, std::string_view b, bool trace = false);

private:
    struct Buffer {
        int* data = nullptr;
        std::size_t bytes = 0;
        void resize(std::size_t count);
        ~Buffer();
    };
    CudaHalo halo_;
    Buffer indices_, table_, boundaries_;
};

} // namespace chps
