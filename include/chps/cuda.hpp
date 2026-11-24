#pragma once

#include "chps/dp.hpp"

namespace chps {

// Reuse exact-sized buffers over serial alignments on the same CUDA device.
class CudaWorkspace {
public:
    CudaWorkspace() = default;
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
    Buffer indices_, table_;
};

} // namespace chps
