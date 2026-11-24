#pragma once
#include "chps/dp.hpp"
#include <array>

namespace chps {
struct SplitStats {
    int proposed_cuts = 0;
    int accepted_cuts = 0;
    int dense_pieces = 0;
    int wavefront_pieces = 0;
    bool global_retry = false;
    bool composition_pruned = false;
};

// Whole-record counts reused across all pairs in a matrix.
class SplitComposition {
public:
    explicit SplitComposition(std::string_view sequence);
    std::size_t lower_bound(const SplitComposition& other) const;
private:
    std::array<std::size_t, 256> counts_{};
};

Result nvpd_split(std::string_view a, std::string_view b, bool trace = false);
Result nvpd_sparse(std::string_view a, std::string_view b, bool trace = false);
enum class SplitPolicy { require_all, retain_verified, single_piece };
Result nvpd_split_profile(std::string_view a, std::string_view b, bool trace, SplitStats& stats,
                          SplitPolicy policy = SplitPolicy::retain_verified,
                          std::size_t lower_bound = 0);
} // namespace chps
