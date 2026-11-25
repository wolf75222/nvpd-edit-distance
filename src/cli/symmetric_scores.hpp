#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chps::cli {

// The diagonal is zero. Each OpenMP task writes only its assigned upper entries.
class SymmetricScores {
public:
    explicit SymmetricScores(std::size_t count) : row_offsets_(count) {
        std::size_t size = 0;
        for (std::size_t row = 0; row < count; ++row) {
            row_offsets_[row] = size;
            const auto width = count - row - 1;
            if (width > scores_.max_size() - size)
                throw std::length_error("distance matrix exceeds addressable storage");
            size += width;
        }
        scores_.resize(size);
    }

    [[nodiscard]] std::size_t pair_count() const noexcept { return scores_.size(); }

    int& upper(std::size_t row, std::size_t column) noexcept {
        return scores_[row_offsets_[row] + (column - row - 1)];
    }

    [[nodiscard]] int score(std::size_t row, std::size_t column) const noexcept {
        if (row == column) return 0;
        const auto [lo, hi] = std::minmax(row, column);
        return scores_[row_offsets_[lo] + (hi - lo - 1)];
    }

private:
    std::vector<std::size_t> row_offsets_;
    std::vector<int> scores_;
};

} // namespace chps::cli
