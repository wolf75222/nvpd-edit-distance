#pragma once

#include "dp_internal.hpp"
#include "chps/fasta.hpp"
#include <algorithm>
#include <numeric>

namespace chps::mpi {

// Each rank shares immutable last-occurrence rows across all of its workers.
class CohortIndexes {
    std::vector<detail::MatchIndex> indexes_;
public:
    explicit CohortIndexes(const std::vector<Record>& records) {
        std::array<bool, 256> seen{};
        std::string alphabet;
        for (const auto& r : records) {
            detail::validate(r.bases, {}, false);
            for (unsigned char c : r.bases) {
                if (!seen[c]) alphabet.push_back(static_cast<char>(c));
                seen[c] = true;
            }
        }
        indexes_.reserve(records.size());
        for (const auto& r : records) indexes_.emplace_back(alphabet, r.bases);
    }

    // Greedily match the shorter string from right to left in the longer one.
    // Keeping the latest occurrence leaves the longest possible prefix for the rest.
    bool contains(std::string_view shorter, std::size_t target, std::size_t length) const {
        const auto& index = indexes_[target];
        const auto stride = length + 1;
        auto end = length;
        for (auto it = shorter.rbegin(); it != shorter.rend(); ++it) {
            const auto symbol = index.symbol[static_cast<unsigned char>(*it)];
            const auto position = index.entries[symbol * stride + end];
            if (!position) return false;
            end = position - 1;
        }
        return true;
    }

    int score(std::string_view a, std::string_view b, std::size_t target,
              std::vector<int>& workspace) const {
        detail::validate(a, b, false);
        const auto stride = b.size() + 1;
        const auto& index = indexes_[target];
        workspace.resize(2 * stride);
        int* previous = workspace.data();
        int* current = previous + stride;
        std::iota(previous, previous + stride, 0);
        int row = 0;
        for (unsigned char c : a) {
            current[0] = ++row;
            const int* last = index.entries.data() + index.symbol[c] * stride;
            #pragma omp simd
            for (int j = 1; j < static_cast<int>(stride); ++j) {
                const int k = last[j];
                current[j] = std::min(1 + std::min(previous[j-1], previous[j]),
                                      previous[std::max(0, k-1)] + (j-k));
            }
            std::swap(previous, current);
        }
        return previous[b.size()];
    }
};
}
