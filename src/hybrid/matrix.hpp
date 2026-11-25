#pragma once

#include "chps/fasta.hpp"

#include <cstddef>
#include <vector>

namespace chps::hybrid {

enum class Policy { edlib, wfa2, gpu, hybrid, hybrid_share, hybrid_counts, cpu_sketch, static_split };

// For a close pair, route to Edlib at or above 5% of the longer sequence.
inline bool prefer_edlib(std::size_t lower_bound, std::size_t longest) {
    return longest != 0 && lower_bound >= longest / 20 + (longest % 20 != 0);
}

struct Matrix {
    std::vector<int> scores;  // Full row-major matrix, including its zero diagonal.
    std::size_t cpu_pairs = 0;
    std::size_t gpu_pairs = 0;
    double sketch_seconds = 0;
    double seconds = 0;
};

Matrix align(const std::vector<Record>& records, Policy policy, int threads);

} // namespace chps::hybrid
