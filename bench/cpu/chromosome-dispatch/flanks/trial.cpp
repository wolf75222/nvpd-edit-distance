#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/split.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdexcept>

namespace chps {
Result original_sparse(std::string_view a, std::string_view b, bool trace);
Result original_split(std::string_view a, std::string_view b, bool trace);
}
namespace {
using Score = chps::Result (*)(std::string_view, std::string_view, bool);
struct Config { const char* name; int threads; Score score; };
}

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::invalid_argument("expected CHROM FASTA1 FASTA2 SCORE");
        const auto left = chps::read_fasta(argv[2]), right = chps::read_fasta(argv[3]);
        if (left.size() != 1 || right.size() != 1) throw std::runtime_error("one record required per FASTA");
        const auto& a = left[0].bases;
        const auto& b = right[0].bases;
        const int expected = std::stoi(argv[4]);
        std::array<Config,8> configs{{
            {"edlib",1,[](auto a, auto b, bool) { return chps::edlib_distance(a,b); }},
            {"wfa2",1,[](auto a, auto b, bool) { return chps::wfa_distance(a,b); }},
            {"split",1,chps::original_split}, {"split",8,chps::original_split},
            {"sparse",1,chps::original_sparse}, {"sparse",8,chps::original_sparse},
            {"reuse",1,chps::nvpd_sparse}, {"reuse",8,chps::nvpd_sparse}}};
        std::mt19937 rng(20261125);
        omp_set_dynamic(0);
        std::cout << "chrom,block,position,method,threads,m,n,score,seconds\n" << std::setprecision(12);
        for (int block = 0; block <= 7; ++block) {
            std::shuffle(configs.begin(),configs.end(),rng);
            int pos = 0;
            for (const auto& c : configs) {
                omp_set_num_threads(c.threads);
                const auto start = std::chrono::steady_clock::now();
                const int score = c.score(a,b,false).distance;
                const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - start;
                if (score != expected) throw std::runtime_error("unexpected chromosome distance");
                std::cout << argv[1] << ',' << block << ',' << ++pos << ',' << c.name << ',' << c.threads
                          << ',' << a.size() << ',' << b.size() << ',' << score << ',' << dt.count() << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
