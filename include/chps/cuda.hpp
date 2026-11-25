#pragma once

#include "chps/dp.hpp"

namespace chps {

// Global uses a uniform left halo. Local handles long reaches per block.
// Automatic also packs short local halos for four-symbol score-only inputs.
enum class CudaHalo { global, local, automatic };
Result nvpd_cuda_local(std::string_view a, std::string_view b, bool trace = false);
// Batch exact score frontiers on cooperatively resident CUDA blocks.
Result nvpd_cuda_persistent(std::string_view a, std::string_view b, bool trace = false);

enum class HybridPairBackend { wfa2, nvpd_cuda, nvpd_cuda_persistent };
constexpr std::string_view hybrid_backend_name(HybridPairBackend backend) {
    if (backend == HybridPairBackend::wfa2) return "wfa2";
    if (backend == HybridPairBackend::nvpd_cuda_persistent) return "nvpd-cuda-persistent";
    return "nvpd-cuda";
}
HybridPairBackend hybrid_pair_backend(std::string_view a, std::string_view b);
Result nvpd_hybrid_score(std::string_view a, std::string_view b, bool trace = false,
                         HybridPairBackend* selected = nullptr);

// Reuse exact-sized buffers over serial alignments on the same CUDA device.
class CudaWorkspace {
public:
    explicit CudaWorkspace(CudaHalo halo = CudaHalo::automatic) : halo_(halo) {}
    CudaWorkspace(const CudaWorkspace&) = delete;
    CudaWorkspace& operator=(const CudaWorkspace&) = delete;
    Result distance(std::string_view a, std::string_view b, bool trace = false);

private:
    friend Result nvpd_cuda_persistent(std::string_view, std::string_view, bool);
    struct Buffer {
        Buffer() = default;
        Buffer(const Buffer&) = delete;
        Buffer& operator=(const Buffer&) = delete;
        ~Buffer() noexcept;

        int* data = nullptr;
        std::size_t bytes = 0;
        void resize(std::size_t count);
    };
    const CudaHalo halo_;
    Buffer indices_, table_, boundaries_;
    std::size_t resident_bytes() const noexcept {
        return indices_.bytes + table_.bytes + boundaries_.bytes;
    }
};

} // namespace chps
