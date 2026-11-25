#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <string_view>

namespace chps::detail::corridor {

// Myers' bit-vector recurrence returns the edit distance to a pattern of at
// most 64 symbols. The bound joins these short alignments at exact anchors.
inline int fragment_score(std::string_view a, std::string_view b) {
    if (a.empty()) return static_cast<int>(b.size());
    std::array<std::uint64_t, 256> mask{};
    for (std::size_t i = 0; i < a.size(); ++i)
        mask[static_cast<unsigned char>(a[i])] |= std::uint64_t(1) << i;
    std::uint64_t pos = ~std::uint64_t(0), neg = 0, hi = std::uint64_t(1) << (a.size() - 1);
    int score = static_cast<int>(a.size());
    for (char c : b) {
        const auto eq = mask[static_cast<unsigned char>(c)], x = eq | neg,
                   d = (((eq & pos) + pos) ^ pos) | eq;
        const auto p = neg | ~(d | pos), m = pos & d;
        score += (p & hi) != 0;
        score -= (m & hi) != 0;
        const auto p2 = (p << 1) | 1, m2 = m << 1;
        pos = m2 | ~(x | p2);
        neg = p2 & x;
    }
    return score;
}

inline int anchored_cost(std::string_view a, std::string_view b, int limit) {
    constexpr std::size_t anchor = 6, radius = 12;
    int cost = 0;
    while (!a.empty() && !b.empty() && cost <= limit) {
        std::size_t same = 0;
        while (same < std::min(a.size(), b.size()) && a[same] == b[same])
            ++same;
        a.remove_prefix(same);
        b.remove_prefix(same);
        if (a.size() <= 64) return cost + fragment_score(a, b);
        if (b.size() < anchor + 1) return cost + fragment_score(b, a);
        auto matches = [&](std::size_t i, std::size_t j) {
            return !std::memcmp(a.data() + i, b.data() + j, anchor);
        };
        if (matches(1, 1)) {
            a.remove_prefix(1);
            b.remove_prefix(1);
            ++cost;
            continue;
        }
        if (matches(0, 1)) {
            b.remove_prefix(1);
            ++cost;
            continue;
        }
        if (matches(1, 0)) {
            a.remove_prefix(1);
            ++cost;
            continue;
        }
        std::size_t aq = 48, bq = std::min(b.size(), aq), consume = 0;
        for (std::size_t candidate : {48u, 56u, 40u, 32u, 24u}) {
            if (candidate + anchor > a.size()) continue;
            for (std::size_t dx = 0; dx <= radius; ++dx) {
                if (candidate + dx + anchor <= b.size() && matches(candidate, candidate + dx)) {
                    aq = candidate;
                    bq = candidate + dx;
                    consume = anchor;
                    break;
                }
                if (dx && candidate - dx + anchor <= b.size() && matches(candidate, candidate - dx)) {
                    aq = candidate;
                    bq = candidate - dx;
                    consume = anchor;
                    break;
                }
            }
            if (consume) break;
        }
        cost += fragment_score(a.substr(0, aq), b.substr(0, bq));
        a.remove_prefix(aq + consume);
        b.remove_prefix(bq + consume);
    }
    return cost > limit ? cost : cost + static_cast<int>(a.size() + b.size());
}

} // namespace chps::detail::corridor
