#pragma once
#include "chps/dp.hpp"

namespace chps {
struct SplitStats {
    int proposed_cuts = 0;
    int accepted_cuts = 0;
    int dense_pieces = 0;
    int wavefront_pieces = 0;
    bool global_retry = false;
};

Result nvpd_split(std::string_view a, std::string_view b, bool trace = false);
enum class SplitPolicy { require_all, retain_verified };
Result nvpd_split_profile(std::string_view a, std::string_view b, bool trace, SplitStats& stats,
                          SplitPolicy policy = SplitPolicy::retain_verified);
} // namespace chps
