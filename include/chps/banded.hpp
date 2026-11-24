#pragma once
#include "chps/dp.hpp"

namespace chps {
Result banded(std::string_view a, std::string_view b, bool trace = false);
Result banded_trim(std::string_view a, std::string_view b, bool trace = false);
}
