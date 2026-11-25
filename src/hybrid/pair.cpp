#include "chps/baselines.hpp"
#include "chps/cuda.hpp"
#include "../cuda/corridor_plan.hpp"
#include <algorithm>
#include <stdexcept>

namespace chps {

HybridPairBackend hybrid_pair_backend(std::string_view a, std::string_view b) {
    constexpr std::size_t max_length = 131072;
    if (std::max(a.size(), b.size()) > max_length) return HybridPairBackend::nvpd_cuda;
    // Canonical orientation also fixes the greedy anchor tie-break on reversals.
    if (a.size() > b.size() || (a.size() == b.size() && a > b)) std::swap(a, b);
    const auto budget = b.size() * 6 / 100;
    if (b.size() - a.size() > budget) return HybridPairBackend::nvpd_cuda;
    const int limit = static_cast<int>(budget);
    return detail::corridor::path_cost(a, b, limit) <= limit
        ? HybridPairBackend::wfa2 : HybridPairBackend::nvpd_cuda;
}

Result nvpd_hybrid_score(std::string_view a, std::string_view b, bool trace,
                         HybridPairBackend* selected) {
    if (trace) throw std::invalid_argument("nvpd-hybrid-score computes the distance only");
    const auto backend = hybrid_pair_backend(a, b);
    if (selected) *selected = backend;
    if (backend == HybridPairBackend::wfa2) return wfa_distance(a, b);
#ifdef CHPS_CUDA
    return nvpd_cuda(a, b);
#else
    throw std::runtime_error("the selected NvPD backend requires a CUDA build");
#endif
}

} // namespace chps
