// Numeric cost in a fixed diagonal band, with an optional optimal-path witness.
#include "edlib.h"

#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <vector>

// The two output integers suffice to certify any band containing this optimal
// path. No aligned strings or operations leave the function.
extern "C" int optimal_witness(const char* a, int m, const char* b, int n, int* reach) {
    const auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
    const auto result = edlibAlign(a, m, b, n, config);
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -2;
    }
    int offset = 0;
    *reach = std::abs(m - n);
    for (int i = 0; i < result.alignmentLength; ++i) {
        offset += (result.alignment[i] == EDLIB_EDOP_INSERT) - (result.alignment[i] == EDLIB_EDOP_DELETE);
        *reach = std::max(*reach, std::abs(offset));
    }
    const int distance = result.editDistance;
    edlibFreeAlignResult(result);
    return distance;
}

extern "C" int fixed_band(const char* a, int m, const char* b, int n, int width) {
    if (width < std::abs(m - n)) return -1;
    const int infinity = m + n + 1;
    const int size = std::min(n + 1, 2 * width + 1);
    std::vector<int> previous(size + 2, infinity), current(size + 2, infinity);
    const int initial_end = std::min(n, width);
    std::iota(previous.begin() + 1, previous.begin() + initial_end + 2, 0);
    int previous_start = 0;
    for (int i = 1; i <= m; ++i) {
        const int start = std::max(0, i - width), end = std::min(n, i + width);
        current.front() = current[end - start + 2] = infinity;
        if (!start) current[1] = i;
        for (int j = std::max(1, start); j <= end; ++j) {
            const int pos = j - start + 1, above = j - previous_start + 1;
            current[pos] = std::min({current[pos - 1] + 1, previous[above] + 1,
                                     previous[above - 1] + (a[i - 1] != b[j - 1])});
        }
        previous.swap(current);
        previous_start = start;
    }
    return previous[n - previous_start + 1];
}
