#include "chps/fasta.hpp"
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace chps {
std::vector<Record> read_fasta(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open FASTA: " + path);
    std::vector<Record> records;
    std::unordered_set<std::string> ids;
    std::string line;
    constexpr std::string_view alphabet = "ACGTRYSWKMBDHVN";
    while (std::getline(file, line)) {
        if (line.starts_with('>')) {
            std::istringstream header(line.substr(1));
            std::string id;
            if (!(header >> id) || !ids.insert(id).second)
                throw std::runtime_error("missing or duplicate FASTA ID: " + id);
            records.push_back({id, {}});
            continue;
        }
        for (unsigned char ch : line) {
            if (std::isspace(ch)) continue;
            char base = static_cast<char>(std::toupper(ch));
            if (records.empty()) throw std::runtime_error("sequence before FASTA header");
            if (alphabet.find(base) == alphabet.npos)
                throw std::runtime_error("invalid DNA symbol in " + records.back().id);
            records.back().bases.push_back(base);
        }
    }
    if (file.bad()) throw std::runtime_error("FASTA read failed: " + path);
    if (records.empty()) throw std::runtime_error("empty FASTA: " + path);
    return records;
}
}
