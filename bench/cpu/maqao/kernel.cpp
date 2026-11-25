#include "chps/dp.hpp"
#include "chps/fasta.hpp"

#include <charconv>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string_view>

namespace {
int positive(std::string_view text) {
    int value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() || value < 1)
        throw std::invalid_argument("expected a positive integer");
    return value;
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 5 && argc != 6)
            throw std::invalid_argument("expected kernel FASTA ROWS THREADS REPEATS [COLUMNS]");
        const auto records = chps::read_fasta(argv[1]);
        const auto rows = static_cast<std::size_t>(positive(argv[2]));
        const auto columns = argc == 6 ? static_cast<std::size_t>(positive(argv[5])) : rows;
        const int threads = positive(argv[3]), repeats = positive(argv[4]);
        if (records.size() != 2 || records[0].bases.size() < rows || records[1].bases.size() < columns)
            throw std::invalid_argument("the two FASTA records must cover ROWS and COLUMNS, respectively");
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
        const auto a = std::string_view(records[0].bases).substr(0, rows);
        const auto b = std::string_view(records[1].bases).substr(0, columns);
        std::cout << "score,seconds\n" << std::setprecision(12);
        for (int i = 0; i < repeats; ++i) {
            const auto start = std::chrono::steady_clock::now();
            const auto score = chps::nvpd_omp(a, b, false).distance;
            const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
            std::cout << score << ',' << elapsed.count() << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
