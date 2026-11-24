#pragma once

#include "chps/dp.hpp"

namespace chps::detail {
// Explicit one-kernel-per-row policy retained for the portable reference measurements.
Result nvpd_kokkos_rows(std::string_view a, std::string_view b, bool trace = false);
Result nvpd_kokkos_opt(std::string_view a, std::string_view b, bool trace = false);
}
