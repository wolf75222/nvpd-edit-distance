#include "chps/baselines.hpp"
#include "dp_internal.hpp"
#include <edlib.h>
#include <cstdlib>
#include <memory>
#include <stdexcept>

extern "C" int chps_wfa_distance(const char*, int, const char*, int);

namespace chps {
Result edlib_distance(std::string_view a, std::string_view b) {
    detail::validate(a, b, false);
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
    auto result = edlibAlign(a.data(), a.size(), b.data(), b.size(), config);
    const int distance = result.status == EDLIB_STATUS_OK ? result.editDistance : -1;
    edlibFreeAlignResult(result);
    if (distance < 0) throw std::runtime_error("edlib alignment failed");
    return {distance, {}, 0};
}
Result edlib_trace(std::string_view a, std::string_view b) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, true);

    // Extended CIGAR describes the query relative to the target. Query B
    // therefore gives insertions and deletions in the project direction A -> B.
    const auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
    const auto result = edlibAlign(b.data(), b.size(), a.data(), a.size(), config);
    const int distance = result.status == EDLIB_STATUS_OK ? result.editDistance : -1;
    char* raw = nullptr;
    if (distance >= 0 && result.alignment)
        raw = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
    edlibFreeAlignResult(result);
    std::unique_ptr<char, decltype(&std::free)> cigar(raw, &std::free);
    if (!cigar) throw std::runtime_error("edlib traceback failed");
    return {distance, cigar.get(), 0};
}

Result wfa_distance(std::string_view a, std::string_view b) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return {static_cast<int>(a.size() + b.size()), {}, 0};
    int distance = chps_wfa_distance(a.data(), a.size(), b.data(), b.size());
    if (distance < 0) throw std::runtime_error("WFA2 alignment did not complete");
    return {distance, {}, 0};
}
}
