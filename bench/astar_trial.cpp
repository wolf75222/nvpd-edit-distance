#include "chps/fasta.hpp"
#include "astarpa.h"

#include <sys/resource.h>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("expected {simple|full} and a FASTA pair");
        const std::string_view mode = argv[1];
        if (mode != "simple" && mode != "full") throw std::invalid_argument("unknown A*PA2 mode");
        const auto records = chps::read_fasta(argv[2]);
        if (records.size() != 2) throw std::invalid_argument("expected two FASTA records");
        for (const auto& r : records)
            if (r.bases.empty() || r.bases.find_first_not_of("ACGT") != std::string::npos)
                throw std::invalid_argument("expected nonempty ACGT sequences");
        const auto& a = records[0].bases;
        const auto& b = records[1].bases;
        auto* align = mode == "simple" ? astarpa2_simple : astarpa2_full;
        std::uint8_t* cigar = nullptr;
        std::uintptr_t cigar_len = 0;
        const auto start = std::chrono::steady_clock::now();
        const auto score = align(reinterpret_cast<const std::uint8_t*>(a.data()), a.size(),
                                 reinterpret_cast<const std::uint8_t*>(b.data()), b.size(),
                                 &cigar, &cigar_len);
        if (!cigar) throw std::runtime_error("missing alignment");
        astarpa_free_cigar(cigar);
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
        usage.ru_maxrss /= 1024;
#endif
        std::cout << std::setprecision(12) << "{\"m\":" << a.size() << ",\"n\":" << b.size()
                  << ",\"distance\":" << score << ",\"seconds\":" << elapsed.count()
                  << ",\"peak_rss_kib\":" << usage.ru_maxrss << "}\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
