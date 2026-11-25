#include "chps/baselines.hpp"
#include "chps/cuda.hpp"
#include "../cuda/corridor_plan.hpp"
#include <algorithm>
#include <stdexcept>

namespace chps {

HybridPairBackend hybrid_pair_backend(std::string_view a, std::string_view b) {
    constexpr std::size_t persistent_min_length = 65536, max_length = 131072;
    if (std::max(a.size(), b.size()) > max_length) return HybridPairBackend::nvpd_cuda;
    // Canonical orientation also fixes the greedy anchor tie-break on reversals.
    if (a.size() > b.size() || (a.size() == b.size() && a > b)) std::swap(a, b);
    const int limit = static_cast<int>(b.size() * 6 / 100);
    const auto difference = b.size() - a.size();
    if (difference <= static_cast<std::size_t>(limit) &&
        (detail::corridor::path_cost(a, b, limit) <= limit ||
         detail::corridor::anchored_cost(a, b, limit) <= limit)) return HybridPairBackend::wfa2;
    if (a.size() < persistent_min_length || difference > b.size() / 4) return HybridPairBackend::nvpd_cuda;
    const auto ambiguous = [](char base) { return detail::corridor::dna_id(base) < 0; };
    if (std::any_of(a.begin(), a.end(), ambiguous) || std::any_of(b.begin(), b.end(), ambiguous))
        return HybridPairBackend::nvpd_cuda;
    return HybridPairBackend::nvpd_cuda_persistent;
}

Result nvpd_hybrid_score(std::string_view a, std::string_view b, bool trace,
                         HybridPairBackend* selected) {
    if (trace) throw std::invalid_argument("nvpd-hybrid-score computes the distance only");
    const auto backend = hybrid_pair_backend(a, b);
    if (selected) *selected = backend;
    if (backend == HybridPairBackend::wfa2) return wfa_distance(a, b);
#ifdef CHPS_CUDA
    if (backend == HybridPairBackend::nvpd_cuda_persistent) return nvpd_cuda_persistent(a, b);
    return nvpd_cuda(a, b);
#else
    throw std::runtime_error("the selected NvPD backend requires a CUDA build");
#endif
}

} // namespace chps
