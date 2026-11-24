#include "chps/baselines.hpp"
#include <edlib.h>
#include <stdexcept>

extern "C" int chps_wfa_distance(const char*, int, const char*, int);

namespace chps {
Result edlib_distance(std::string_view a, std::string_view b) {
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
    auto result = edlibAlign(a.data(), a.size(), b.data(), b.size(), config);
    int distance = result.editDistance;
    bool ok = result.status == EDLIB_STATUS_OK && distance >= 0;
    edlibFreeAlignResult(result);
    if (!ok) throw std::runtime_error("edlib alignment failed");
    return {distance, {}, 0};
}
Result wfa_distance(std::string_view a, std::string_view b) {
    if (a.empty() || b.empty()) return {static_cast<int>(a.size() + b.size()), {}, 0};
    int distance = chps_wfa_distance(a.data(), a.size(), b.data(), b.size());
    if (distance < 0) throw std::runtime_error("WFA2 alignment did not complete");
    return {distance, {}, 0};
}
}
