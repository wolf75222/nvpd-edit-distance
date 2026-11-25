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
namespace {
struct EdlibAlignment {
    explicit EdlibAlignment(EdlibAlignResult result) : value(result) {}
    ~EdlibAlignment() { edlibFreeAlignResult(value); }
    EdlibAlignment(const EdlibAlignment&) = delete;
    EdlibAlignment& operator=(const EdlibAlignment&) = delete;
    const EdlibAlignResult value;
};

Result wfa(std::string_view a, std::string_view b, bool trace, bool bidirectional) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    // The upstream traceback buffer is indexed by 2*(m+n) in a signed int.
    const std::size_t combined = a.size() + b.size();
    if (trace && combined > INT_MAX / 2)
        throw std::length_error("WFA2 traceback requires combined length <= INT_MAX/2");
    // WFA2 2.3.6 grows its signed sequence buffer by half after adding three 64-byte pads.
    const std::size_t padded = combined + 3 * 64;
    if (padded + padded / 2 > INT_MAX)
        throw std::length_error("WFA2 sequence buffer exceeds its signed index range");
    char* raw = nullptr;
    const int distance = chps_wfa_align(a.data(), static_cast<int>(a.size()),
                                        b.data(), static_cast<int>(b.size()),
                                        bidirectional, trace ? &raw : nullptr);
    std::unique_ptr<char, decltype(&std::free)> cigar(raw, &std::free);
    if (distance < 0) throw std::runtime_error("WFA2 alignment did not complete");
    Result result{distance, cigar ? cigar.get() : "", 0};
    std::replace(result.cigar.begin(), result.cigar.end(), 'M', '=');
    return result;
}
} // namespace

Result edlib_distance(std::string_view a, std::string_view b) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, false);
    const auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
    int distance;
    {
        const EdlibAlignment alignment(edlibAlign(a.data(), static_cast<int>(a.size()),
                                                  b.data(), static_cast<int>(b.size()), config));
        distance = alignment.value.status == EDLIB_STATUS_OK ? alignment.value.editDistance : -1;
    }
    if (distance < 0) throw std::runtime_error("edlib alignment failed");
    return {distance, {}, 0};
}

Result edlib_trace(std::string_view a, std::string_view b) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, true);

    // Edlib's query must be B for the extended CIGAR to describe A -> B.
    const auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
    int distance;
    std::unique_ptr<char, decltype(&std::free)> cigar(nullptr, &std::free);
    {
        const EdlibAlignment alignment(edlibAlign(b.data(), static_cast<int>(b.size()),
                                                  a.data(), static_cast<int>(a.size()), config));
        const auto& result = alignment.value;
        if (result.status != EDLIB_STATUS_OK || result.editDistance < 0)
            throw std::runtime_error("edlib alignment failed");
        distance = result.editDistance;
        if (!result.alignment) throw std::runtime_error("edlib traceback failed");
        cigar.reset(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED));
    }
    if (!cigar) throw std::runtime_error("edlib traceback failed");
    return {distance, cigar.get(), 0};
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
