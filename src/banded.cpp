#include "chps/banded.hpp"
#include <algorithm>
#include <climits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace chps {
namespace {
// All paths of cost <= k stay inside |i-j| <= k. The two sentinels
// allow the left edge of the band to move by one column between rows.
int bounded(std::string_view a, std::string_view b, int k) {
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const int width = std::min(n + 1, 2 * k + 1), inf = k + 1;
    std::vector<int> previous(width + 2), current(width + 2);
    int prev_lo = 0, prev_hi = std::min(n, k);
    std::iota(previous.begin() + 1, previous.begin() + prev_hi + 2, 0);
    previous[0] = previous[prev_hi + 2] = inf;
    for (int i = 1; i <= m; ++i) {
        const int lo = std::max(0, i - k), hi = std::min(n, i + k);
        current[0] = current[hi - lo + 2] = inf;
        if (lo == 0) current[1] = i;
        int best = lo == 0 ? i : inf;
        for (int j = std::max(1, lo); j <= hi; ++j) {
            const int x = j - lo + 1, y = j - prev_lo + 1;
            const int v = std::min({current[x - 1] + 1, previous[y] + 1,
                                    previous[y - 1] + (a[i - 1] != b[j - 1])});
            current[x] = v;
            best = std::min(best, v);
        }
        if (best > k) return inf;
        previous.swap(current);
        prev_lo = lo;
    }
    return previous[n - prev_lo + 1];
}
}

Result banded(std::string_view a, std::string_view b, bool trace) {
    if (trace) throw std::invalid_argument("banded traceback unavailable");
    if (a.size() < b.size()) std::swap(a, b);
    if (a.size() > (INT_MAX - 3) / 2) throw std::length_error("sequence exceeds band index range");
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    if (!n) return {m, {}, 0};
    for (int k = std::max(1, m - n); ; k = std::min(m, 2 * k)) {
        const int score = bounded(a, b, k);
        if (score <= k) return {score, {}, 0};
    }
}

Result banded_trim(std::string_view a, std::string_view b, bool trace) {
    std::size_t prefix = 0;
    while (prefix < std::min(a.size(), b.size()) && a[prefix] == b[prefix]) ++prefix;
    a.remove_prefix(prefix);
    b.remove_prefix(prefix);
    while (!a.empty() && !b.empty() && a.back() == b.back()) {
        a.remove_suffix(1);
        b.remove_suffix(1);
    }
    return banded(a, b, trace);
}
}
