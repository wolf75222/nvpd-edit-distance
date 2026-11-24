#include "chps/packed.hpp"
#include "chps/baselines.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <limits>
#include <stdexcept>
#if defined(__x86_64__) && defined(__BMI2__)
#include <immintrin.h>
#endif

namespace chps {
namespace {
constexpr auto code = [] {
    std::array<unsigned char, 256> table{};
    table.fill(4);
    table['A'] = 0;
    table['C'] = 1;
    table['G'] = 3;
    table['T'] = 2;
    return table;
}();

int byte_extension(std::string_view a, std::string_view b, int i, int j, int length) {
    int offset = 0;
    while (length - offset >= 8) {
        std::uint64_t x, y;
        std::memcpy(&x, a.data() + i + offset, sizeof x);
        std::memcpy(&y, b.data() + j + offset, sizeof y);
        const auto difference = x ^ y;
        if (difference) {
            if constexpr (std::endian::native == std::endian::little)
                return offset + std::countr_zero(difference) / 8;
            if constexpr (std::endian::native == std::endian::big)
                return offset + std::countl_zero(difference) / 8;
            break;
        }
        offset += 8;
    }
    while (offset < length && a[i + offset] == b[j + offset]) ++offset;
    return offset;
}

int reverse_bytes(std::string_view a, std::string_view b, int end_a, int end_b, int length) {
    int offset = 0;
    while (length - offset >= 8) {
        std::uint64_t x, y;
        std::memcpy(&x, a.data() + end_a - offset - 8, sizeof x);
        std::memcpy(&y, b.data() + end_b - offset - 8, sizeof y);
        const auto difference = x ^ y;
        if (difference) {
            if constexpr (std::endian::native == std::endian::little)
                return offset + std::countl_zero(difference) / 8;
            if constexpr (std::endian::native == std::endian::big)
                return offset + std::countr_zero(difference) / 8;
            break;
        }
        offset += 8;
    }
    while (offset < length && a[end_a-1-offset] == b[end_b-1-offset]) ++offset;
    return offset;
}

template<bool Guarded, class Extend>
int frontier(std::string_view a, std::string_view b, Extend extend) {
    constexpr int limit = packed_front_limit;
    if (std::max(a.size(), b.size()) > static_cast<std::size_t>(std::numeric_limits<int>::max() - limit - 1))
        throw std::length_error("input exceeds signed alignment coordinate range");
    const auto n = static_cast<int>(a.size()), m = static_cast<int>(b.size());
    if (!n || !m) return std::max(n, m);
    if (std::abs(n - m) > limit) return biwfa(a, b).distance;
    std::array<int, 2 * limit + 3> front_a, front_b;
    int* prev = front_a.data() + limit + 1;
    int* next = front_b.data() + limit + 1;
    prev[0] = extend(0, 0, std::min(n, m));
    if (n == m && prev[0] == n) return 0;
    int low = 0, high = 0;
    for (int d = 1; d <= limit; ++d) {
        const int new_low = Guarded ? std::max(-m, -d) : -d;
        const int new_high = Guarded ? std::min(n, d) : d;
        if constexpr (!Guarded) {
            prev[-d-1] = prev[-d] = prev[d] = prev[d+1] = -2;
        }
        for (int k = new_low; k <= new_high; ++k) {
            int x = -1;
            if constexpr (Guarded) {
                if (k >= low && k <= high && prev[k] >= 0 && prev[k] < n && prev[k] - k < m)
                    x = prev[k] + 1;
                if (k > low && k - 1 <= high && prev[k-1] >= 0 && prev[k-1] < n)
                    x = std::max(x, prev[k-1] + 1);
                if (k < high && k + 1 >= low && prev[k+1] >= 0 && prev[k+1] - (k+1) < m)
                    x = std::max(x, prev[k+1]);
            } else {
                x = std::max({prev[k]+1, prev[k-1]+1, prev[k+1]});
            }
            const int remaining = std::min(n-x, m-(x-k));
            if (x >= 0 && remaining > 0) x += extend(x, x-k, remaining);
            next[k] = x;
        }
        if (n-m >= new_low && n-m <= new_high && next[n-m] >= n) return d;
        std::swap(prev, next);
        low = new_low;
        high = new_high;
    }
    return biwfa(a, b).distance;
}

template<class Forward, class Reverse>
int bidirectional(std::string_view a, std::string_view b, Forward forward, Reverse reverse) {
    constexpr int limit = packed_front_limit;
    if (std::max(a.size(), b.size()) > static_cast<std::size_t>(std::numeric_limits<int>::max() - limit - 1))
        throw std::length_error("input exceeds signed alignment coordinate range");
    const auto n = static_cast<int>(a.size()), m = static_cast<int>(b.size());
    if (!n || !m) return std::max(n, m);
    const int delta = n - m;
    if (std::abs(delta) > limit) return biwfa(a, b).distance;
    constexpr int padding = limit / 2 + 2;
    using Layer = std::array<int, 2 * padding + 1>;
    Layer first_f, first_r, scratch;
    int* f = first_f.data()+padding;
    int* r = first_r.data()+padding;
    int* next = scratch.data()+padding;
    const int common = std::min(n,m);
    f[0] = forward(0, 0, common);
    if (delta == 0 && f[0] == n) return 0;
    r[0] = reverse(n, m, common);
    auto step = [&](auto&& extend, int d, int*& current, const int* other, int other_d) {
        current[-d-1] = current[-d] = current[d] = current[d+1] = -2;
        const int intersect_lo = std::max(-d, delta-other_d);
        const int intersect_hi = std::min(d, delta+other_d);
        for (int k = -d; k <= d; ++k) {
            int x = std::max({current[k]+1, current[k-1]+1, current[k+1]});
            const int room = std::min(n-x, m-(x-k));
            if (room > 0) x += extend(x,x-k,room);
            next[k] = x;
            if (k >= intersect_lo && k <= intersect_hi && x >= n-other[delta-k]) return true;
        }
        std::swap(current,next);
        return false;
    };
    auto from_end = [&](int i, int j, int len) { return reverse(n-i, m-j, len); };
    for (int d = 1; d <= limit/2; ++d) {
        if (step(forward,d,f,r,d-1)) return 2*d-1;
        if (step(from_end,d,r,f,d)) return 2*d;
    }
    return biwfa(a,b).distance;
}

}

PackedSequence::PackedSequence(std::string_view bases) : bases_(bases) {
    unsigned invalid = 0;
    #pragma omp simd reduction(|:invalid)
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const char c = bases[i];
        invalid |= static_cast<unsigned>(!((c == 'A') | (c == 'C') | (c == 'G') | (c == 'T')));
    }
    if (invalid) return;
    words_.resize((bases.size() + 31) / 32 + 1);
    std::size_t i = 0;
#if defined(__x86_64__) && defined(__BMI2__)
    for (; bases.size() - i >= 32; i += 32) {
        std::uint64_t out = 0;
        for (int j = 0; j < 4; ++j) {
            std::uint64_t lane;
            std::memcpy(&lane, bases.data()+i+8*j, sizeof lane);
            out |= _pext_u64(lane, UINT64_C(0x0606060606060606)) << (16*j);
        }
        words_[i/32] = out;
    }
#endif
    for (; i < bases.size(); ++i) {
        const auto c = code[static_cast<unsigned char>(bases[i])];
        words_[i / 32] |= static_cast<std::uint64_t>(c) << (2 * (i % 32));
    }
}

std::uint64_t PackedSequence::block(std::size_t pos) const {
    const auto index = pos / 32, shift = 2 * (pos % 32);
    const auto lo = words_[index] >> shift;
    return shift ? lo | (words_[index+1] << (64-shift)) : lo;
}

int packed_score(const PackedSequence& a, const PackedSequence& b, PackedFront front) {
    if (!a.is_acgt() || !b.is_acgt()) return biwfa(a.bases(), b.bases()).distance;
    auto forward = [&](int i, int j, int length) {
        int offset = byte_extension(a.bases(), b.bases(), i, j, std::min(8,length));
        if (offset != 8) return offset;
        for (; offset < length; offset += 32) {
            const auto difference = a.block(i + offset) ^ b.block(j + offset);
            if (difference) return std::min(length, offset + std::countr_zero(difference) / 2);
        }
        return length;
    };
    auto reverse = [&](int i, int j, int length) {
        int offset = reverse_bytes(a.bases(), b.bases(), i, j, std::min(8,length));
        if (offset != 8) return offset;
        while (length-offset >= 32) {
            const auto difference = a.block(i-offset-32) ^ b.block(j-offset-32);
            if (difference) return offset + std::countl_zero(difference) / 2;
            offset += 32;
        }
        return offset + reverse_bytes(a.bases(), b.bases(), i-offset, j-offset, length-offset);
    };
    if (front == PackedFront::guarded) return frontier<true>(a.bases(), b.bases(), forward);
    if (front == PackedFront::one_sided) return frontier<false>(a.bases(), b.bases(), forward);
    return bidirectional(a.bases(), b.bases(), forward, reverse);
}

int byte_wavefront_score(std::string_view a, std::string_view b) {
    return bidirectional(a, b, [&](int i, int j, int len) { return byte_extension(a, b, i, j, len); },
                         [&](int i, int j, int len) { return reverse_bytes(a, b, i, j, len); });
}

Result packed_wavefront(std::string_view a, std::string_view b, bool trace) {
    if (trace) return biwfa(a, b, true);
    return {packed_score(PackedSequence(a), PackedSequence(b)), {}, 0};
}
}
