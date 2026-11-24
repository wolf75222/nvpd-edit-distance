#include "chps/adaptive.hpp"
#include "dp_internal.hpp"
#include "nvpd_band_row.hpp"
#include <array>
#include <climits>
#include <numeric>
#include <stdexcept>

namespace chps {
namespace {
int bounded(std::string_view a, std::string_view b, const detail::MatchIndex& matches,
            int threshold, AdaptiveWork& work) {
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const int delta = n - m, margin = (threshold - delta) / 2;
    const int right = delta + margin, width = right + margin + 1;
    const int inf = threshold + 1;
    std::vector<int> prev(width + 2), curr(width + 2);
    int prev_lo = 0, prev_hi = right;
    std::iota(prev.begin() + 1, prev.begin() + prev_hi + 2, 0);
    for (int i = 1; i <= m; ++i) {
        const int lo = std::max(0, i - margin), hi = std::min(n, i + right);
        prev[0] = prev[prev_hi-prev_lo+2] = inf;
        if (!lo) curr[1] = i;
        const auto symbol = matches.symbol[static_cast<unsigned char>(a[i - 1])];
        const int* last = matches.entries.data() + symbol * (b.size() + 1);
        const int best = detail::nvpd_band_row(prev.data(), curr.data(), last, lo, hi, prev_lo, threshold);
        work.cells += hi - std::max(1, lo) + 1;
        if (best > threshold && (lo || i > threshold)) return inf;
        prev.swap(curr);
        prev_lo = lo;
        prev_hi = hi;
    }
    return prev[n-prev_lo+1];
}
}

int nvpd_adaptive_score(std::string_view a, std::string_view b, AdaptiveWork& work) {
    work = {};
    detail::validate(a, b, false);
    if (std::max(a.size(), b.size()) > (INT_MAX - 3) / 2)
        throw std::length_error("sequences exceed band index range");
    if (a.size() > b.size()) std::swap(a, b);
    while (static_cast<std::size_t>(work.prefix) < a.size() && a[work.prefix] == b[work.prefix])
        ++work.prefix;
    a.remove_prefix(work.prefix);
    b.remove_prefix(work.prefix);
    while (!a.empty() && a.back() == b.back()) {
        a.remove_suffix(1);
        b.remove_suffix(1);
        ++work.suffix;
    }
    if (a.empty()) return static_cast<int>(b.size());
    // Literal symbol deficits bound the number of required edits from below.
    // Substituting aligned mismatches and inserting the remainder gives a path.
    std::array<int, 256> balance{};
    int upper = static_cast<int>(b.size() - a.size());
    for (std::size_t j = 0; j < a.size(); ++j) {
        --balance[static_cast<unsigned char>(a[j])];
        ++balance[static_cast<unsigned char>(b[j])];
        upper += a[j] != b[j];
    }
    for (std::size_t j = a.size(); j < b.size(); ++j) ++balance[static_cast<unsigned char>(b[j])];
    int lower = 0;
    for (int count : balance) lower += std::max(0, count);
    if (lower == upper) return lower;
    const detail::MatchIndex matches(a, b);
    for (int k = std::max(1, lower); ; k = std::min(upper, 2*k)) {
        ++work.passes;
        work.threshold = k;
        const int score = bounded(a, b, matches, k, work);
        if (score <= k) return score;
    }
}

Result nvpd_adaptive(std::string_view a, std::string_view b, bool trace) {
    if (trace) throw std::invalid_argument("nvpd-adaptive supports score only");
    AdaptiveWork work;
    return {nvpd_adaptive_score(a, b, work), {}, 0};
}
}
