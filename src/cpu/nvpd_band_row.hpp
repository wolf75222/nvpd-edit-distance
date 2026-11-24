#pragma once
#include <algorithm>

namespace chps::detail {
// Index 0 of prev is the sentinel immediately left of the previous interval.
inline int nvpd_band_row(const int* prev, int* curr, const int* last, int lo,
                         int hi, int prev_lo, int threshold) {
    const int inf = threshold + 1;
    int best = inf;
    #pragma omp simd reduction(min:best)
    for (int j = std::max(1, lo); j <= hi; ++j) {
        const int k = last[j];
        const int p = j - prev_lo;
        const int via_match = prev[std::max(0, k - prev_lo)] + (j - k);
        const int value = std::min(inf, std::min(1 + std::min(prev[p], prev[p+1]), via_match));
        curr[j-lo+1] = value;
        best = std::min(best, value);
    }
    return best;
}
}
