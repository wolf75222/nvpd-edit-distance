#pragma once

#include "chps/dp.hpp"

#include <array>
#include <vector>

namespace chps::detail {

void validate(std::string_view a, std::string_view b, bool trace);
Result empty_result(std::string_view a, std::string_view b, bool trace);
std::string traceback(const int* table, std::string_view a, std::string_view b);

// Sadiq et al., Eq. (2). k is the most recent match in B[1..j].
// Matching at k costs prev[k-1], followed by j-k insertions.
#ifdef __CUDACC__
__host__ __device__
#endif
inline int nvpd_cell(const int* prev, int j, int k) {
    if (k == j) return prev[j-1];
    const int upper = 1 + (prev[j] < prev[j-1] ? prev[j] : prev[j-1]);
    if (k == 0) return upper;
    const int via_match = prev[k-1] + (j - k);
    return upper < via_match ? upper : via_match;
}

// A row per distinct symbol of A. Each entry gives its last 1-based
// position in the corresponding prefix of B, or zero when it is absent.
struct MatchIndex {
    MatchIndex(std::string_view a, std::string_view b);

    std::array<int, 256> symbol{};
    std::vector<int> entries;
};

} // namespace chps::detail
