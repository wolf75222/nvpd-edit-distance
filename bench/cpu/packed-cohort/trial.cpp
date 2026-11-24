#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/packed.hpp"
#include <omp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace {
using Clock = std::chrono::steady_clock;
using Matrix = std::vector<int>;
using chps::Record;

std::size_t pair_index(std::size_t n, std::size_t i, std::size_t j) {
    return i * (2 * n-i-1) / 2 + j-i-1;
}

class Trial {
public:
    const std::vector<Record>& input;
    std::vector<std::string_view> rows;
    std::vector<std::size_t> row_id;
    std::vector<chps::PackedSequence> packed;

    Trial(const std::vector<Record>& input, bool collapse, bool encode) : input(input) {
        std::unordered_map<std::string_view, std::size_t> index;
        for (const auto& record : input) {
            if (collapse) {
                const auto [it, fresh] = index.emplace(record.bases, rows.size());
                row_id.push_back(it->second);
                if (!fresh) continue;
            }
            rows.push_back(record.bases);
        }
        if (encode) {
            packed.resize(rows.size());
            #pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < rows.size(); ++i) packed[i] = chps::PackedSequence(rows[i]);
        }
    }

    Matrix compute(const std::string& method) const {
        const auto n = rows.size();
        Matrix scores(n * (n-1) / 2);
        #pragma omp parallel for schedule(dynamic, 1)
        for (std::size_t i=0; i<n; ++i) {
            for (std::size_t j=i+1; j<n; ++j) {
                int distance;
                if (method == "packed") distance = chps::packed_score(packed[i], packed[j]);
                else if (method == "guarded" || method == "one-sided")
                    distance = chps::packed_score(packed[i], packed[j], method == "guarded"
                                                  ? chps::PackedFront::guarded : chps::PackedFront::one_sided);
                else if (method == "byte") distance = chps::byte_wavefront_score(rows[i], rows[j]);
                else if (method == "wfa2") distance = chps::wfa_distance(rows[i], rows[j]).distance;
                else distance = chps::edlib_distance(rows[i], rows[j]).distance;
                scores[pair_index(n,i,j)] = distance;
            }
        }
        if (row_id.empty()) return scores;
        const auto m = input.size();
        Matrix output(m * (m-1) / 2);
        #pragma omp parallel for schedule(static)
        for (std::size_t i=0; i<m; ++i) {
            for (std::size_t j=i+1; j<m; ++j) {
                const auto x = row_id[i], y = row_id[j];
                output[pair_index(m,i,j)] = x == y ? 0 : scores[pair_index(n,std::min(x,y),std::max(x,y))];
            }
        }
        return output;
    }
};

void check() {
    std::mt19937 rng(1001);
    const std::string alphabet = "ACGT";
    std::vector<std::string> short_dna = {""};
    for (std::size_t i=0; short_dna[i].size()<5; ++i)
        for (char c : std::string("AC")) short_dna.push_back(short_dna[i]+c);
    auto check_pair = [](const std::string& a, const std::string& b) {
        const auto oracle = chps::edlib_distance(a,b).distance;
        const chps::PackedSequence x(a), y(b);
        for (auto front : {chps::PackedFront::guarded, chps::PackedFront::one_sided, chps::PackedFront::bidirectional})
            if (chps::packed_score(x,y,front) != oracle) throw std::runtime_error("incorrect packed front");
        if (chps::byte_wavefront_score(a,b) != oracle) throw std::runtime_error("incorrect byte front");
    };
    for (const auto& a : short_dna) for (const auto& b : short_dna) check_pair(a,b);
    for (std::size_t k=0; k<1024; ++k) {
        std::string a(128 + k%32, 'A');
        for (auto& c : a) c = alphabet[rng()%4];
        auto b = a.substr(k%32);
        b.insert(0, (k/32)%32, 'T');
        for (std::size_t z=0; z<k%19; ++z) b[rng()%b.size()] = alphabet[rng()%4];
        check_pair(a,b);
        check_pair(b,a);
    }
    for (int edits = 88; edits <= 104; ++edits) {
        std::string a(193, 'A'), b = a;
        std::vector<int> sites(a.size());
        std::iota(sites.begin(),sites.end(),0);
        std::shuffle(sites.begin(),sites.end(),rng);
        for (int i = 0; i < edits; ++i) b[sites[i]] = alphabet[1+rng()%3];
        check_pair(a,b);
        std::string random = a;
        for (auto& c : random) c = alphabet[rng()%4];
        check_pair(b,random);
        check_pair(random.substr(0,random.size()-edits),random);
    }
    check_pair(std::string(200,'A'),std::string(200,'T'));
    check_pair(std::string(300,'C'),"C");
    check_pair("ACGTRYSWKMBDHVN","ACGTNYSWRMBDHVN");
    check_pair("ACGTN","ACGTA");
}
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
    if (argc != 6) throw std::invalid_argument("trial FASTA OUTPUT_DIR THREADS REPETITIONS MODE[full|scale|weak]");
    const std::string dir = argv[2];
    const int threads = std::stoi(argv[3]), reps=std::stoi(argv[4]);
    const std::string mode = argv[5];
    if (mode != "full" && mode != "scale" && mode != "weak") throw std::invalid_argument("unknown mode");
    const bool scale = mode != "full";
    auto input = chps::read_fasta(argv[1]);
    if (mode == "weak") {
        const double target = static_cast<double>(input.size()) * (input.size()-1) * threads / 64;
        input.resize(static_cast<std::size_t>(std::lround((1 + std::sqrt(1 + 4*target))/2)));
    }
    omp_set_num_threads(threads);
    if (input.size()<2) throw std::runtime_error("expected at least two records");
    Matrix oracle;
    const char* cache = std::getenv("CHPS_ORACLE_CACHE");
    if (scale && cache) {
        std::ifstream file(cache, std::ios::binary);
        std::size_t n;
        if (!file.read(reinterpret_cast<char*>(&n), sizeof n) || n < input.size())
            throw std::runtime_error("invalid Edlib oracle cache");
        Matrix full(n*(n-1)/2);
        if (!file.read(reinterpret_cast<char*>(full.data()), full.size()*sizeof(int)))
            throw std::runtime_error("truncated Edlib oracle cache");
        oracle.reserve(input.size()*(input.size()-1)/2);
        for (std::size_t i=0; i<input.size(); ++i)
            for (std::size_t j=i+1; j<input.size(); ++j) oracle.push_back(full[pair_index(n,i,j)]);
    } else {
        const Trial truth(input,true,false);
        oracle=truth.compute("edlib");
        if (!scale && cache) {
            std::ofstream file(cache, std::ios::binary);
            const auto n = input.size();
            file.write(reinterpret_cast<const char*>(&n), sizeof n);
            file.write(reinterpret_cast<const char*>(oracle.data()), oracle.size()*sizeof(int));
            if (!file) throw std::runtime_error("failed to save the Edlib oracle");
        }
    }
    if (!scale) {
        std::map<int,std::uint64_t> frequency;
        for (int value : oracle) ++frequency[value];
        std::ofstream hist(dir+"/distances.csv");
        hist << "distance,pairs\n";
        for (auto [d,count] : frequency) hist << d << ',' << count << '\n';
    }
    std::vector<std::pair<std::string,bool>> cases;
    for (const auto* method : {"edlib","wfa2","byte","packed","guarded","one-sided"}) {
        if (scale && std::string(method) != "packed") continue;
        cases.emplace_back(method,false);
        if (mode != "weak") cases.emplace_back(method,true);
    }
    std::mt19937 rng(1001);
    std::cout << "trial,threads,method,collapse,input_rows,unique_rows,kernel_calls,encoding_bytes,prepare_s,compute_s,total_s\n";
    for (int rep=0; rep<=reps; ++rep) {
        std::shuffle(cases.begin(),cases.end(),rng);
        for (const auto& [method,collapse] : cases) {
            const auto before = Clock::now();
            Trial trial(input,collapse,method == "packed" || method == "guarded" || method == "one-sided");
            const auto ready = Clock::now();
            auto got = trial.compute(method);
            const auto done = Clock::now();
            if (got != oracle) throw std::runtime_error("matrix disagrees with Edlib");
            const auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
            std::size_t bytes=0;
            for (const auto& row : trial.packed) bytes+=row.encoding_bytes();
            const auto n=trial.rows.size();
            std::cout << rep << ',' << threads << ',' << method << ',' << collapse << ','
                      << input.size() << ',' << n << ',' << n*(n-1)/2 << ',' << bytes << ','
                      << seconds(before,ready) << ',' << seconds(ready,done) << ',' << seconds(before,done) << std::endl;
        }
    }
}
