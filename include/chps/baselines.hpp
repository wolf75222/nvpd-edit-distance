#pragma once
#include "chps/dp.hpp"

namespace chps {
Result edlib_distance(std::string_view a, std::string_view b);
Result edlib_trace(std::string_view a, std::string_view b);
Result wfa_distance(std::string_view a, std::string_view b);
Result wfa_trace(std::string_view a, std::string_view b);
Result biwfa(std::string_view a, std::string_view b, bool trace = false);
}
