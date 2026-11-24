#include "sparse_frontier.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <stdexcept>

namespace chps::detail {
namespace {
constexpr int max_budget = 128;

int extend(std::string_view a, std::string_view b, int i, int j) {
    const int available = std::min(int(a.size()) - i, int(b.size()) - j);
    int n = 0;
    if (available && a[i] == b[j]) {
        for (int block : {4096, 64})
            while (available - n >= block && std::memcmp(a.data()+i+n, b.data()+j+n, block) == 0)
                n += block;
        while (n < available && a[i+n] == b[j+n]) ++n;
    }
    return i + n;
}
}

int sparse_frontier(std::string_view a, std::string_view b, int budget) {
    if (budget < 0 || budget > max_budget)
        throw std::invalid_argument("sparse frontier budget must lie in [0,128]");
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    if (std::abs(m - n) > budget) return -1;
    std::array<int, 2*max_budget+3> row0, row1;
    row0.fill(-1);
    row1.fill(-1);
    int *prev = row0.data()+max_budget+1, *next = row1.data()+max_budget+1;
    prev[0] = extend(a, b, 0, 0);
    if (m == n && prev[0] == m) return 0;
    // On diagonal k = i-j, keep the furthest query offset after extending exact matches.
    for (int s = 1; s <= budget; ++s) {
        for (int k = -s; k <= s; ++k) {
            int i = -1;
            const int del = prev[k-1], sub = prev[k], ins = prev[k+1];
            if (del >= 0 && del < m) i = del + 1;
            if (sub >= 0 && sub < m && sub-k < n) i = std::max(i, sub+1);
            if (ins >= 0 && ins-(k+1) < n) i = std::max(i, ins);
            next[k] = i < 0 ? -1 : extend(a, b, i, i-k);
        }
        if (std::abs(m-n) <= s && next[m-n] == m) return s;
        std::swap(prev, next);
    }
    return -1;
}
}
