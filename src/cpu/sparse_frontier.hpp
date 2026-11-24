#pragma once

#include <string_view>

namespace chps::detail {
// Return -1 if the unit-cost distance exceeds the supplied budget.
int sparse_frontier(std::string_view a, std::string_view b, int budget);
}
