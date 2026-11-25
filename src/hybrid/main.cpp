#include "matrix.hpp"
#include "sketch.hpp"

#include <charconv>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace {
using chps::Record;
using chps::hybrid::Matrix;
using chps::hybrid::Policy;
constexpr auto usage =
    "editdist-route --input FILE --policy cpu-edlib|cpu-wfa2|gpu|hybrid|hybrid-share|cpu-sketch|static "
    "[--threads N] [--warmups N] [--repeats N] [--matrix FILE] [--decisions FILE]";

int number(std::string_view value, int minimum) {
    int n = 0;
    auto [end, err] = std::from_chars(value.data(), value.data() + value.size(), n);
    if (err != std::errc{} || end != value.data() + value.size() || n < minimum)
        throw std::invalid_argument("invalid count: " + std::string(value));
    return n;
}

bool same_file(const std::string& a, const std::string& b) {
    namespace fs = std::filesystem;
    if (a.empty() || b.empty()) return false;
    if (fs::exists(a) && fs::exists(b)) return fs::equivalent(a, b);
    return fs::weakly_canonical(fs::absolute(a)) == fs::weakly_canonical(fs::absolute(b));
}

std::ofstream output(const std::string& path) {
    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.open(path);
    return out;
}

void write_matrix(const std::string& path, const std::vector<Record>& records, const Matrix& data) {
    auto out = output(path);
    out << "id";
    for (const auto& record : records) out << '\t' << record.id;
    out << '\n';
    const auto n = records.size();
    for (std::size_t i = 0; i < n; ++i) {
        out << records[i].id;
        for (std::size_t j = 0; j < n; ++j) out << '\t' << data.scores[i * n + j];
        out << '\n';
    }
    out.close();
}

void write_decisions(const std::string& path, const std::vector<Record>& records, const Matrix& data) {
    std::vector<chps::hybrid::Sketch> words;
    words.reserve(records.size());
    for (const auto& r : records) words.push_back(chps::hybrid::sketch(r.bases));
    auto out = output(path);
    out << "i\tj\tid_i\tid_j\tshared\tsampled\tjaccard\tclass\tdistance\n" << std::setprecision(12);
    for (std::size_t i = 0; i < records.size(); ++i) {
        for (std::size_t j = i + 1; j < records.size(); ++j) {
            const auto match = chps::hybrid::overlap(words[i], words[j]);
            out << i << '\t' << j << '\t' << records[i].id << '\t' << records[j].id << '\t'
                << match.common << '\t' << match.sampled << '\t';
            if (match.sampled == 0) out << "NA\tfallback";
            else out << static_cast<double>(match.common) / match.sampled << '\t' << (match.close() ? "close" : "far");
            out << '\t' << data.scores[i * records.size() + j] << '\n';
        }
    }
    out.close();
}
}

int main(int argc, char** argv) {
    try {
        const std::unordered_map<std::string, Policy> policies{
            {"gpu", Policy::gpu}, {"cpu-edlib", Policy::edlib}, {"cpu-wfa2", Policy::wfa2},
            {"hybrid", Policy::hybrid}, {"cpu-sketch", Policy::cpu_sketch},
            {"static", Policy::static_split}, {"hybrid-share", Policy::hybrid_share}};
        int threads = 1, warmups = 0, repeats = 1;
        std::string input, name, matrix_path, decisions_path;
        for (int i = 1; i < argc; ++i) {
            const std::string_view key(argv[i]);
            if (key == "--help") { std::cout << usage << '\n'; return 0; }
            if (++i == argc) throw std::invalid_argument(usage);
            const std::string_view value(argv[i]);
            if (key == "--input") input = value;
            else if (key == "--policy") name = value;
            else if (key == "--threads") threads = number(value, 1);
            else if (key == "--warmups") warmups = number(value, 0);
            else if (key == "--repeats") repeats = number(value, 1);
            else if (key == "--matrix") matrix_path = value;
            else if (key == "--decisions") decisions_path = value;
            else throw std::invalid_argument(usage);
        }
        if (input.empty() || !policies.contains(name)) throw std::invalid_argument(usage);
        if (same_file(matrix_path, input) || same_file(decisions_path, input) || same_file(matrix_path, decisions_path))
            throw std::invalid_argument("input and output paths must be distinct");
        omp_set_dynamic(0);
        const auto records = chps::read_fasta(input);
        for (const auto& r : records)
            if (r.bases.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::length_error("sequence exceeds the alignment index range");
        std::cout.exceptions(std::ios::badbit | std::ios::failbit);
        std::cout << "phase,iteration,policy,threads,sequences,pairs,cpu_pairs,gpu_pairs,sketch_seconds,seconds\n"
                  << std::setprecision(12);
        const auto policy = policies.at(name);
        for (bool warmup : {true, false}) {
            const int count = warmup ? warmups : repeats;
            for (int it = 0; it < count; ++it) {
                const auto result = chps::hybrid::align(records, policy, threads);
                std::cout << (warmup ? "warmup" : "measure") << ',' << it + 1 << ',' << name << ',' << threads
                          << ',' << records.size() << ',' << result.cpu_pairs + result.gpu_pairs << ','
                          << result.cpu_pairs << ',' << result.gpu_pairs << ',' << result.sketch_seconds << ','
                          << result.seconds << '\n' << std::flush;
                if (!warmup && it + 1 == count) {
                    if (!matrix_path.empty()) write_matrix(matrix_path, records, result);
                    if (!decisions_path.empty()) write_decisions(decisions_path, records, result);
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "editdist-route: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
