#pragma once
#include "chps/dp.hpp"
#include <charconv>
#include <stdexcept>

inline void replay(std::string_view a, std::string_view b, const chps::Result& result) {
    std::string_view cigar = result.cigar;
    size_t i = 0, j = 0;
    int cost = 0;
    while (!cigar.empty()) {
        int count = 0;
        const auto [op, ec] = std::from_chars(cigar.data(), cigar.data() + cigar.size(), count);
        if (ec != std::errc{} || op == cigar.data() + cigar.size() || count <= 0)
            throw std::runtime_error("malformed CIGAR");
        cigar.remove_prefix(op + 1 - cigar.data());
        for (int k = 0; k < count; ++k) {
            switch (*op) {
            case '=': if (a.at(i++) != b.at(j++)) throw std::runtime_error("bad match"); break;
            case 'X': if (a.at(i++) == b.at(j++)) throw std::runtime_error("bad mismatch"); ++cost; break;
            case 'I': (void)b.at(j++); ++cost; break;
            case 'D': (void)a.at(i++); ++cost; break;
            default: throw std::runtime_error("bad CIGAR operation");
            }
        }
    }
    if (i != a.size() || j != b.size() || cost != result.distance)
        throw std::runtime_error("traceback failed to consume sequences or achieve optimal cost");
}

