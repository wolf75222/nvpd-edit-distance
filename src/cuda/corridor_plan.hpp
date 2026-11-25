#pragma once

#include "corridor_anchor.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <string_view>
#include <vector>

namespace chps::detail::corridor {
#ifndef CHPS_NVPD_CORRIDOR_ROWS
#define CHPS_NVPD_CORRIDOR_ROWS 16
#endif
#ifndef CHPS_NVPD_CORRIDOR_COLUMNS
#define CHPS_NVPD_CORRIDOR_COLUMNS 128
#endif
#ifndef CHPS_NVPD_CORRIDOR_ANCHORS
#define CHPS_NVPD_CORRIDOR_ANCHORS 1
#endif
constexpr int rows = CHPS_NVPD_CORRIDOR_ROWS;
constexpr int columns = CHPS_NVPD_CORRIDOR_COLUMNS;
constexpr int threads = std::min(columns,128);
static_assert(rows >= 1 && rows <= 64);
static_assert(columns >= 32 && columns <= 512 && columns % 32 == 0);

inline int dna_id(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Construct a feasible alignment. A short exact anchor chooses an indel;
// a failed anchor search costs one substitution. The cost is an upper bound.
inline int path_cost(std::string_view a, std::string_view b, int limit) {
    constexpr std::size_t anchor = 6, lookahead = 8;
    std::size_t i = 0, j = 0;
    int cost = 0;
    while (i < a.size() && j < b.size() && cost <= limit) {
        if (a[i] == b[j]) { ++i; ++j; continue; }
        bool advanced = false;
        if (a.size() - i >= anchor && b.size() - j >= anchor) {
            auto match = [&](std::size_t x, std::size_t y) {
                return std::memcmp(a.data() + x, b.data() + y, anchor) == 0;
            };
            if (a.size() - i > anchor && b.size() - j > anchor && match(i+1,j+1)) {
                ++i; ++j; ++cost;
                continue;
            }
            for (std::size_t s = 1; s <= lookahead; ++s) {
                if (s + anchor <= b.size() - j && match(i,j+s)) {
                    j += s; cost += static_cast<int>(s); advanced = true; break;
                }
                if (s + anchor <= a.size() - i && match(i+s,j)) {
                    i += s; cost += static_cast<int>(s); advanced = true; break;
                }
            }
        }
        if (!advanced) { ++i; ++j; ++cost; }
    }
    return cost > limit ? cost : cost + static_cast<int>(std::max(a.size()-i,b.size()-j));
}

struct Chunk {
    unsigned mask[4]{};
    int before[4]{};
};
static_assert(sizeof(Chunk) == 32);

#ifdef __CUDACC__
__host__ __device__
#endif
inline int last(const Chunk& c, int j, int symbol) {
    const unsigned mask = c.mask[symbol] & ((2u << (j & 31)) - 1u);
    if (!mask) return c.before[symbol];
#ifdef __CUDA_ARCH__
    return (j & ~31) + 31 - __clz(mask);
#else
    return (j & ~31) + static_cast<int>(std::bit_width(mask)) - 1;
#endif
}

inline std::vector<Chunk> chunk_index(std::string_view b) {
    std::vector<Chunk> chunks(b.size()/32 + 1);
    std::array<int,4> latest{};
    for (int j = 1; j <= static_cast<int>(b.size()); ++j) {
        const int id = dna_id(b[j-1]);
        if (id < 0) return {};
        auto& chunk = chunks[j/32];
        if ((j & 31) == 0) std::copy(latest.begin(),latest.end(),chunk.before);
        chunk.mask[id] |= 1u << (j & 31);
        latest[id] = j;
    }
    return chunks;
}

enum class Path { size, upper_bound, length_difference, alphabet, wide_halo, tiles };
struct Layout {
    Path path = Path::size;
    int upper = -1, pitch = 0, capacity = 0;
    std::vector<Chunk> index;
    std::vector<int> edges;
};

// A sweep can compile the CUDA translation unit with a different tile geometry.
static inline Layout make_layout(std::string_view a, std::string_view b) {
    Layout layout;
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    if (m < 32768 || m > n || n-m > n/4) return layout;
    layout.upper = path_cost(a,b,n/4);
    if constexpr (CHPS_NVPD_CORRIDOR_ANCHORS) {
        if (layout.upper != n-m)
            layout.upper = std::min(layout.upper,anchored_cost(a,b,std::min(layout.upper,n/4)));
    }
    if (layout.upper > n/4) { layout.path = Path::upper_bound; return layout; }
    if (layout.upper == n-m) { layout.path = Path::length_difference; return layout; }
    layout.path = Path::alphabet;
    if (std::any_of(a.begin(),a.end(),[](char c){return dna_id(c)<0;})) return layout;
    layout.index = chunk_index(b);
    if (layout.index.empty()) return layout;
    const int blocks = n/columns+1;
    layout.edges.resize(static_cast<std::size_t>(blocks)*rows);
    for (int block = 0; block < blocks; ++block) {
        int* edge = layout.edges.data() + static_cast<std::size_t>(block)*rows;
        const int begin = block*columns, end = std::min(n,begin+columns-1);
        edge[rows-1] = begin;
        for (int t = rows-1; t > 0; --t) {
            const int j = edge[t];
            int p = std::max(0,j-1);
            for (int s = 0; s < 4; ++s) {
                const int k = last(layout.index[j/32],j,s);
                if (k) p = std::min(p,k-1);
            }
            edge[t-1] = p;
        }
        layout.pitch = std::max(layout.pitch,end-edge[0]+1);
        // Exceptional symbol gaps can make redundant halo work dominate.
        if (layout.pitch > 8*columns) { layout.path = Path::wide_halo; return layout; }
    }
    layout.path = Path::tiles;
    layout.capacity = std::min(n+1,layout.upper+1+2*(columns-1));
    return layout;
}

struct Slice {
    int first, last;
};

static inline Slice output_slice(int end_row, int m, int n, int upper) {
    const int margin = (upper-(n-m))/2, right = n-m+margin;
    const int first = (std::max(0,end_row-margin)/columns)*columns;
    const int last_block = std::min(n,end_row+right)/columns;
    return {first,std::min(n,(last_block+1)*columns-1)};
}
} // namespace chps::detail::corridor
