#pragma once
#include "chps/dp.hpp"
#include <cstdint>

namespace chps {
struct AdaptiveWork {
    std::uint64_t cells = 0;
    int passes = 0, threshold = 0, prefix = 0, suffix = 0;
};
Result nvpd_adaptive(std::string_view a, std::string_view b, bool trace = false);
int nvpd_adaptive_score(std::string_view a, std::string_view b, AdaptiveWork& work);
}
