#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace chps::hybrid {

using Sketch = std::vector<std::uint64_t>;
constexpr std::size_t sketch_size = 256;
constexpr unsigned word_length = 15;

// SplitMix64 finalizer. Its odd multipliers and xor shifts permute the word codes.
inline std::uint64_t word_hash(std::uint64_t x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

inline Sketch sketch(std::string_view seq) {
    static constexpr auto codes = [] {
        std::array<unsigned, 256> values{};
        values.fill(4);
        values['A'] = 0; values['C'] = 1; values['G'] = 2; values['T'] = 3;
        return values;
    }();
    constexpr auto mask = (UINT64_C(1) << (2 * word_length)) - 1;
    Sketch heap;
    heap.reserve(sketch_size);
    std::unordered_set<std::uint64_t> present;
    present.reserve(sketch_size + 1);
    std::uint64_t word = 0;
    unsigned filled = 0;
    for (unsigned char c : seq) {
        const auto code = codes[c];
        if (code == 4) {
            filled = 0;
            word = 0;
            continue;
        }
        word = ((word << 2) | code) & mask;
        if (filled < word_length && ++filled < word_length) continue;
        const auto hash = word_hash(word);
        if (heap.size() == sketch_size && hash >= heap.front()) continue;
        if (!present.insert(hash).second) continue;
        if (heap.size() == sketch_size) {
            std::pop_heap(heap.begin(), heap.end());
            present.erase(heap.back());
            heap.pop_back();
        }
        heap.push_back(hash);
        std::push_heap(heap.begin(), heap.end());
    }
    std::sort_heap(heap.begin(), heap.end());
    return heap;
}

// Among the bottom hashes of the union, count those present in both sets.
struct Overlap {
    std::size_t common = 0, sampled = 0;
    bool close() const { return sampled > 0 && 4 * common >= sampled; }
};

inline Overlap overlap(const Sketch& a, const Sketch& b) {
    if (a.empty() || b.empty()) return {};
    std::size_t i = 0, j = 0, sampled = 0, common = 0;
    while (sampled < sketch_size && (i < a.size() || j < b.size())) {
        if (j == b.size() || (i < a.size() && a[i] < b[j])) ++i;
        else if (i == a.size() || b[j] < a[i]) ++j;
        else {
            ++i;
            ++j;
            ++common;
        }
        ++sampled;
    }
    return {common, sampled};
}

} // namespace chps::hybrid
