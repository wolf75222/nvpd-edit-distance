#include "chps/fasta.hpp"
#include <algorithm>
#include <array>
#include <cctype>
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

bool is_upper_acgt(std::string_view line) {
    unsigned other = 0;
    #pragma omp simd reduction(|:other)
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        other |= (c != 'A') & (c != 'C') & (c != 'G') & (c != 'T');
    }
    return other == 0;
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
