#include "chps/fasta.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <filesystem>
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

// FAI offsets refer to uncompressed bytes. Samtools can omit empty records.
class FastaIndex {
    std::ifstream file_;
    std::uintmax_t file_size_;
    std::size_t remaining_ = 0, line_bases_ = 0, line_bytes_ = 0;
    std::string pending_;

    static std::size_t number(std::string_view field) {
        std::size_t n = 0;
        const auto [end, error] = std::from_chars(field.data(), field.data() + field.size(), n);
        if (field.empty() || error != std::errc{} || end != field.data() + field.size())
            throw std::runtime_error("invalid FAI integer");
        return n;
    }

public:
    FastaIndex(const std::string& fasta, const std::string& index) : file_(index) {
        if (!file_) throw std::runtime_error("cannot open FAI: " + index);
        file_size_ = std::filesystem::file_size(fasta);
    }

    std::size_t next(std::string_view id, std::size_t offset) {
        if (remaining_) throw std::runtime_error("FAI length exceeds decoded sequence");
        std::string row;
        row.swap(pending_);
        if (row.empty() && !std::getline(file_, row)) {
            if (!file_.eof() || file_.bad()) throw std::runtime_error("FAI read failed");
            return 0;
        }
        if (!row.empty() && row.back() == '\r') row.pop_back();
        std::string_view rest(row);
        std::array<std::string_view, 5> fields;
        for (std::size_t i = 0; i < fields.size() - 1; ++i) {
            const auto tab = rest.find('\t');
            if (tab == rest.npos) throw std::runtime_error("FAI requires five tab-separated fields");
            fields[i] = rest.substr(0, tab);
            rest.remove_prefix(tab + 1);
        }
        fields.back() = rest;
        if (fields[0] != id) {
            // Keep this entry for the next header. line() will reject any bases
            // in the current record because its permitted length remains zero.
            pending_.swap(row);
            return 0;
        }
        const auto length = number(fields[1]);
        if (number(fields[2]) != offset) throw std::runtime_error("FAI offset disagrees with FASTA");
        line_bases_ = number(fields[3]);
        line_bytes_ = number(fields[4]);
        if (offset > file_size_ || length > file_size_ - offset || line_bases_ > length ||
            (length && !line_bases_) || line_bytes_ < line_bases_ || line_bytes_ - line_bases_ > 2 ||
            (!length && line_bytes_))
            throw std::runtime_error("invalid FAI length or line geometry");
        remaining_ = length;
        return length;
    }

    void line(std::size_t decoded, std::size_t physical, std::size_t terminator) {
        if (!remaining_ && !decoded && physical == terminator) return;
        const auto expected = std::min(remaining_, line_bases_);
        if (decoded != expected || physical - terminator != decoded ||
            (remaining_ > line_bases_ && physical != line_bytes_) ||
            (remaining_ == line_bases_ && terminator && physical != line_bytes_))
            throw std::runtime_error("FAI length or line geometry disagrees with FASTA");
        remaining_ -= decoded;
    }

    void finish() {
        if (remaining_ || !pending_.empty() || file_.peek() != std::char_traits<char>::eof() || file_.bad())
            throw std::runtime_error("FAI has incorrect length or extra entries");
    }
};

template<bool Indexed>
std::vector<Record> decode(const std::string& path, FastaIndex* index) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open FASTA: " + path);
    std::vector<Record> records;
    std::unordered_set<std::string> ids;
    std::string line;
    std::size_t offset = 0;
    while (std::getline(file, line)) {
        std::size_t bytes = 0, terminator = 0;
        if constexpr (Indexed) {
            bytes = line.size() + !file.eof();
            terminator = !file.eof() + static_cast<std::size_t>(line.ends_with('\r'));
            offset += bytes;
        }
        if (line.starts_with('>')) {
            std::istringstream header(line.substr(1));
            std::string id;
            if (!(header >> id) || !ids.insert(id).second)
                throw std::runtime_error("missing or duplicate FASTA ID: " + id);
            records.push_back({id, {}});
            if constexpr (Indexed) records.back().bases.reserve(index->next(id, offset));
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
        if constexpr (Indexed) index->line(line.size(), bytes, terminator);
        records.back().bases.append(line);
    }
    if (file.bad()) throw std::runtime_error("FASTA read failed: " + path);
    if (records.empty()) throw std::runtime_error("empty FASTA: " + path);
    if constexpr (Indexed) index->finish();
    return records;
}
} // namespace

std::vector<Record> read_fasta(const std::string& path) {
    return decode<false>(path, nullptr);
}

std::vector<Record> read_fasta(const std::string& path, const std::string& index_path) {
    FastaIndex index(path, index_path);
    return decode<true>(path, &index);
}
}
