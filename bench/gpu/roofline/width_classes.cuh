#pragma once

#include "local_tiles.cuh"
#include <algorithm>
#include <vector>

namespace chps::detail {
struct WidthClasses {
    int small_count = 0, small_pitch = 0, large_pitch = 0;
};

// Isolate the widest tenth of blocks, retaining the original output position
// in the last boundary. All blocks with the same pitch stay in the same class.
inline WidthClasses group_widths(LocalTilePlan &plan, int n, int columns) {
    const int blocks = static_cast<int>(plan.left.size() / local_tile_rows);
    if (!blocks)
        return {};
    std::vector<int> widths(blocks);
    for (int block = 0; block < blocks; ++block) {
        const auto *edge = plan.left.data() + block * local_tile_rows;
        int start = edge[local_tile_rows - 1];
        widths[block] = start - edge[0] + std::min(columns, n - start + 1);
    }
    auto sorted = widths;
    auto cutoff = sorted.begin() + (9LL * blocks - 1) / 10;
    std::nth_element(sorted.begin(), cutoff, sorted.end());
    WidthClasses g;
    std::vector<int> ordered(plan.left.size());
    int a = 0, b = static_cast<int>(
                   std::count_if(widths.begin(), widths.end(), [&](int w) { return w <= *cutoff; }));
    g.small_count = b;
    for (int block = 0; block < blocks; ++block) {
        int where;
        if (widths[block] <= *cutoff) {
            where = a++;
            g.small_pitch = std::max(g.small_pitch, widths[block]);
        } else {
            where = b++;
            g.large_pitch = std::max(g.large_pitch, widths[block]);
        }
        std::copy_n(plan.left.data() + block * local_tile_rows, local_tile_rows,
                    ordered.data() + where * local_tile_rows);
    }
    plan.left.swap(ordered);
    return g;
}
} // namespace chps::detail
