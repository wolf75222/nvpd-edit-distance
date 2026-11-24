#include "dp_internal.hpp"

#include <algorithm>
#include <climits>
#include <numeric>
#include <cstdint>
#include <omp.h>
#include <type_traits>
#include <stdexcept>

namespace chps::detail {

void validate(std::string_view a, std::string_view b, bool trace) {
    if (a.size() >= INT_MAX || b.size() >= INT_MAX || a.size() + b.size() >= INT_MAX)
        throw std::length_error("combined length must be less than INT_MAX");
    if (trace && std::max(a.size(), b.size()) > 20'000)
        throw std::length_error("traceback supports at most 20000 bases per input");
}

Result empty_result(std::string_view a, std::string_view b, bool trace) {
    const int cost = static_cast<int>(std::max(a.size(), b.size()));
    std::string cigar;
    if (trace && cost > 0) cigar = std::to_string(cost) + (a.empty() ? 'I' : 'D');
    return {cost, cigar, 0};
}

std::string traceback(const int* table, std::string_view a, std::string_view b) {
    const std::size_t stride = b.size() + 1;
    std::size_t i = a.size(), j = b.size();
    std::string ops;
    while (i || j) {
        const auto here = i * stride + j;
        if (i && j && table[here] == table[here - stride - 1] + (a[i-1] != b[j-1])) {
            ops += a[i-1] == b[j-1] ? '=' : 'X';
            --i;
            --j;
        } else if (i && table[here] == table[here - stride] + 1) {
            ops += 'D';
            --i;
        } else {
            ops += 'I';
            --j;
        }
    }
    std::string cigar;
    for (std::size_t end = ops.size(); end > 0;) {
        const char op = ops[end - 1];
        std::size_t start = end - 1;
        while (start > 0 && ops[start - 1] == op) --start;
        cigar += std::to_string(end - start);
        cigar += op;
        end = start;
    }
    return cigar;
}

MatchIndex::MatchIndex(std::string_view a, std::string_view b) {
    symbol.fill(-1);
    std::vector<unsigned char> alphabet;
    for (unsigned char c : a) {
        if (symbol[c] >= 0) continue;
        symbol[c] = static_cast<int>(alphabet.size());
        alphabet.push_back(c);
    }
    const std::size_t stride = b.size() + 1;
    entries.resize(alphabet.size() * stride);
    #pragma omp parallel for schedule(static)
    for (int s = 0; s < static_cast<int>(alphabet.size()); ++s) {
        int* row = entries.data() + s * stride;
        int last = 0;
        for (int j = 1; j <= static_cast<int>(b.size()); ++j) {
            if (static_cast<unsigned char>(b[j - 1]) == alphabet[s]) last = j;
            row[j] = last;
        }
    }
}

} // namespace chps::detail

namespace chps {

Result sequential(std::string_view a, std::string_view b, bool trace) {
    detail::validate(a, b, trace);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    if (!trace && a.size() < b.size()) std::swap(a, b);
    const std::size_t stride = b.size() + 1;
    std::vector<int> table((trace ? a.size() + 1 : 2) * stride);
    std::iota(table.begin(), table.begin() + stride, 0);
    int* prev = table.data();
    for (int i = 1; i <= static_cast<int>(a.size()); ++i) {
        int* curr = table.data() + (trace ? i : i % 2) * stride;
        curr[0] = i;
        for (int j = 1; j <= static_cast<int>(b.size()); ++j) {
            curr[j] = std::min({prev[j] + 1, curr[j-1] + 1,
                                prev[j-1] + (a[i-1] != b[j-1])});
        }
        prev = curr;
    }
    return {prev[b.size()], trace ? detail::traceback(table.data(), a, b) : "", 0};
}

Result antidiagonal(std::string_view a, std::string_view b) {
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, false);
    if (a.size() > b.size()) std::swap(a, b);
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const std::size_t stride = a.size() + 1;
    std::vector<int> table(3 * stride);
    #pragma omp parallel
    {
        for (int d = 0; d <= m + n; ++d) {
            int* curr = table.data() + (d % 3) * stride;
            const int* prev = table.data() + ((d % 3 + 2) % 3) * stride;
            const int* diag = table.data() + ((d % 3 + 1) % 3) * stride;
            const int first = std::max(0, d - n), last = std::min(m, d);
            #pragma omp for schedule(static)
            for (int i = first; i <= last; ++i) {
                const int j = d - i;
                if (i == 0 || j == 0) curr[i] = d;
                else curr[i] = std::min({prev[i-1] + 1, prev[i] + 1,
                                         diag[i-1] + (a[i-1] != b[j-1])});
            } // All producers finish before the next diagonal reuses a buffer.
        }
    }
    return {table[((m + n) % 3) * stride + m], "", 0};
}

Result nvpd_omp(std::string_view a, std::string_view b, bool trace) {
    detail::validate(a, b, trace);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    if (!trace && a.size() > b.size()) std::swap(a, b);
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const std::size_t stride = b.size() + 1;
    detail::MatchIndex matches(a, b);
    std::vector<std::uint8_t> offsets(matches.entries.size());
    bool compact = true;
    for (std::size_t s = 0; compact && s < matches.entries.size(); s += stride) {
        for (int j = 1; j <= n; ++j) {
            int gap = j - matches.entries[s + j];
            if (gap > UINT8_MAX) { compact = false; break; }
            offsets[s + j] = static_cast<std::uint8_t>(gap);
        }
    }
    if (compact) std::vector<int>().swap(matches.entries);
    else std::vector<std::uint8_t>().swap(offsets);
    std::vector<int> table((trace ? a.size() + 1 : 2) * stride);
    std::iota(table.begin(), table.begin() + stride, 0);
    const auto fill = [&](const auto& entries) {
        using Code = typename std::decay_t<decltype(entries)>::value_type;
        #pragma omp parallel
        {
            for (int i = 1; i <= m; ++i) {
                const int* prev = table.data() + (trace ? i - 1 : (i - 1) % 2) * stride;
                int* curr = table.data() + (trace ? i : i % 2) * stride;
                const auto s = matches.symbol[static_cast<unsigned char>(a[i-1])];
                const auto* last = entries.data() + s * stride;
                if (omp_get_thread_num() == 0) curr[0] = i;
                #pragma omp for simd schedule(static)
                for (int j = 1; j <= n; ++j) {
                    int k = last[j];
                    if constexpr (std::is_same_v<Code, std::uint8_t>) k = j - k;
                    // k=j gives the matching diagonal. If k=0, this candidate
                    // is (i-1)+j >= max(i,j), a valid upper bound on the cost.
                    const int via_match = prev[std::max(0, k - 1)] + (j - k);
                    curr[j] = std::min(1 + std::min(prev[j-1], prev[j]), via_match);
                } // One barrier per row; the next row reads only this row.
            }
        }
    };
    if (compact) fill(offsets);
    else fill(matches.entries);
    return {table[(trace ? m : m % 2) * stride + n],
            trace ? detail::traceback(table.data(), a, b) : "", 0};
}

} // namespace chps
