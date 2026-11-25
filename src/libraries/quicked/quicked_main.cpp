#include "chps/fasta.hpp"
#include "chps/quicked.h"
#include <sys/resource.h>
#include <chrono>
#include <climits>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("usage: editdist-quicked PAIR.fa");
        const auto records = chps::read_fasta(argv[1]);
        if (records.size() != 2) throw std::invalid_argument("expected two records");
        const auto& a = records[0].bases;
        const auto& b = records[1].bases;
        for (const auto& r : records) {
            if (r.bases.find_first_not_of("ACGT") != std::string::npos)
                throw std::invalid_argument("QuickEd comparison requires ACGT");
            if (r.bases.size() > INT_MAX / 100)
                throw std::length_error("QuickEd percentage arithmetic exceeds safe range");
        }
        const auto start = std::chrono::steady_clock::now();
        const auto result = chps_quicked_certified(a.data(), static_cast<int>(a.size()), b.data(), static_cast<int>(b.size()));
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        if (result.distance < 0) throw std::runtime_error("QuickEd failed to return a distance");
        rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
        usage.ru_maxrss /= 1024;
#endif
        const char* verifier = a.empty() || b.empty() ? "length" : "edlib";
        std::cout << std::setprecision(12) << "{\"algorithm\":\"quicked\",\"verified_by\":\"" << verifier << "\",\"m\":" << a.size()
                  << ",\"method\":\"quicked-certified\",\"source\":\"" << result.source
                  << "\",\"first_score\":" << result.first_score << ",\"retries\":" << result.retries
                  << ",\"band_percent\":" << result.band_percent
                  << ",\"n\":" << b.size() << ",\"distance\":" << result.distance << ",\"seconds\":" << elapsed.count()
                  << ",\"peak_rss_kib\":" << usage.ru_maxrss << "}\n";
        if (!std::cout.flush()) throw std::runtime_error("cannot write result to standard output");
    } catch (const std::exception& e) {
        std::cerr << "editdist-quicked: " << e.what() << '\n';
        return 1;
    }
}
