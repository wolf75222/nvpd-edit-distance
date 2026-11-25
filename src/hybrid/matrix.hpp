#pragma once

#include "chps/fasta.hpp"

#include <cstddef>
#include <vector>

namespace chps::hybrid {

enum class Policy { edlib, wfa2, gpu, hybrid, hybrid_share, cpu_sketch, static_split };

struct Matrix {
    std::vector<int> scores;  // Full row-major matrix, including its zero diagonal.
    std::size_t cpu_pairs = 0;
    std::size_t gpu_pairs = 0;
    double sketch_seconds = 0;
    double seconds = 0;
};

Matrix align(const std::vector<Record>& records, Policy policy, int threads);

} // namespace chps::hybrid
