#include "chps/baselines.hpp"
#include "dp_internal.hpp"
#include <edlib.h>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <memory>
#include <stdexcept>

extern "C" int chps_wfa_align(const char*, int, const char*, int, bool, char**);

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

static Result wfa(std::string_view a, std::string_view b, bool trace, bool bidirectional) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    // The upstream traceback buffer is indexed by 2*(m+n) in a signed int.
    if (trace && a.size() + b.size() > INT_MAX / 2)
        throw std::length_error("WFA2 traceback requires combined length <= INT_MAX/2");
    char* raw = nullptr;
    const int distance = chps_wfa_align(a.data(), a.size(), b.data(), b.size(),
                                        bidirectional, trace ? &raw : nullptr);
    std::unique_ptr<char, decltype(&std::free)> cigar(raw, &std::free);
    if (distance < 0) throw std::runtime_error("WFA2 alignment did not complete");
    Result result{distance, cigar ? cigar.get() : "", 0};
    std::replace(result.cigar.begin(), result.cigar.end(), 'M', '=');
    return result;
}
Result wfa_distance(std::string_view a, std::string_view b) {
    return wfa(a, b, false, false);
}
Result wfa_trace(std::string_view a, std::string_view b) {
    return wfa(a, b, true, false);
}
Result biwfa(std::string_view a, std::string_view b, bool trace) {
    return wfa(a, b, trace, true);
}
}
