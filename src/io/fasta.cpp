#include "chps/fasta.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace chps {
namespace {
constexpr auto base_codes = [] {
    std::array<char, 256> codes{};
    for (unsigned char base : std::string_view("ACGTRYSWKMBDHVN")) {
        codes[base] = static_cast<char>(base);
        codes[base + ('a' - 'A')] = static_cast<char>(base);
    }
    return codes;
}();

#if defined(__GNUC__) || defined(__clang__)
bool acgt_block(const char* text) {
    using Bytes = unsigned char __attribute__((vector_size(16)));
    Bytes bases;
    std::memcpy(&bases, text, sizeof(bases));
    const Bytes mask = (bases == 'A') | (bases == 'C') | (bases == 'G') | (bases == 'T');
    std::array<std::uint64_t, 2> halves;
    std::memcpy(halves.data(), &mask, sizeof(mask));
    return (halves[0] & halves[1]) == UINT64_MAX;
}
#endif

bool is_upper_acgt(std::string_view line) {
#if defined(__GNUC__) || defined(__clang__)
    if (line.size() >= 16) {
        const char* last = line.data() + line.size() - 16;
        for (const char* block = line.data(); block < last; block += 16)
            if (!acgt_block(block)) return false;
        // The final block may overlap its predecessor and stays inside the line.
        return acgt_block(last);
    }
#endif
    return std::all_of(line.begin(), line.end(), [](char c) {
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    });
}

} // namespace

std::vector<Record> read_fasta(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open FASTA: " + path);
    std::vector<Record> records;
    std::unordered_set<std::string> ids;
    std::string line;
    while (std::getline(file, line)) {
        if (line.starts_with('>')) {
            std::istringstream header(line.substr(1));
            std::string id;
            if (!(header >> id) || !ids.insert(id).second)
                throw std::runtime_error("missing or duplicate FASTA ID: " + id);
            records.push_back({id, {}});
            continue;
        }
        if (records.empty()) {
            if (!std::all_of(line.begin(), line.end(), [](unsigned char c) { return std::isspace(c); }))
                throw std::runtime_error("sequence before FASTA header");
            continue;
        }
        if (!is_upper_acgt(line)) {
            std::size_t used = 0;
            for (unsigned char ch : line) {
                if (const char base = base_codes[ch])
                    line[used++] = base;
                else if (!std::isspace(ch))
                    throw std::runtime_error("invalid DNA symbol in " + records.back().id);
            }
            line.resize(used);
        }
        records.back().bases.append(line);
    }
    if (file.bad()) throw std::runtime_error("FASTA read failed: " + path);
    if (records.empty()) throw std::runtime_error("empty FASTA: " + path);
    return records;
}
}
