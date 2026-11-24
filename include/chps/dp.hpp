#pragma once

#include <cstddef>
#include <string>
#include <string_view>

namespace chps {

struct Result {
    int distance;
    std::string cigar;             // Forward A -> B, with runs of =, X, I, D.
    std::size_t device_bytes;      // Peak explicit CUDA allocations, in bytes.
};

Result sequential(std::string_view a, std::string_view b, bool trace = false);
Result antidiagonal(std::string_view a, std::string_view b);
Result nvpd_omp(std::string_view a, std::string_view b, bool trace = false);
Result nvpd_cuda(std::string_view a, std::string_view b, bool trace = false);

} // namespace chps
