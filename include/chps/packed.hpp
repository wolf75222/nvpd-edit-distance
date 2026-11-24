#pragma once

#include "chps/dp.hpp"
#include <cstdint>
#include <vector>

namespace chps {

// The original string must outlive this view. An empty word array marks IUPAC input.
class PackedSequence {
public:
    PackedSequence() = default;
    explicit PackedSequence(std::string_view bases);
    std::string_view bases() const { return bases_; }
    bool is_acgt() const { return !words_.empty(); }
    std::size_t encoding_bytes() const { return words_.size() * sizeof(std::uint64_t); }
    std::uint64_t block(std::size_t pos) const;

private:
    std::string_view bases_;
    std::vector<std::uint64_t> words_;
};

inline constexpr int packed_front_limit = 96;
enum class PackedFront { bidirectional, one_sided, guarded };
int packed_score(const PackedSequence& a, const PackedSequence& b, PackedFront front = PackedFront::bidirectional);
int byte_wavefront_score(std::string_view a, std::string_view b);
Result packed_wavefront(std::string_view a, std::string_view b, bool trace = false);

}
