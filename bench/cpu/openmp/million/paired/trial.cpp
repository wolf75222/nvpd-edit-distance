#include "chps/dp.hpp"
#include "chps/fasta.hpp"

#include <array>
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::invalid_argument("expected FASTA PAIR SCORE");
        const auto records = chps::read_fasta(argv[1]);
        const int pair = std::stoi(argv[2]), expected = std::stoi(argv[3]);
        if (records.size() != 2 || pair < 1 || pair > 5) throw std::invalid_argument("unexpected input");
        const auto& a = records[0].bases;
        const auto& b = records[1].bases;
        constexpr std::array methods{"diagonal", "nvpd-omp"};
        std::cout << "pair,block,position,method,length_a,length_b,score,seconds\n" << std::setprecision(12);
        for (int block = 0; block <= 3; ++block) {
            for (int pos = 0; pos < 2; ++pos) {
                const int index = (block + pair + pos) % 2;
                const auto start = std::chrono::steady_clock::now();
                const auto score = index ? chps::nvpd_omp(a, b).distance : chps::antidiagonal(a, b).distance;
                const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
                if (score != expected) throw std::runtime_error("disagreement with Edlib/WFA2");
                std::cout << pair << ',' << block << ',' << pos+1 << ',' << methods[index] << ','
                          << a.size() << ',' << b.size() << ',' << score << ',' << elapsed.count() << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
