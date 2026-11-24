#include "chps/fasta.hpp"
#include <edlib.h>
#include <sys/resource.h>
#include <chrono>
#include <climits>
#include <iomanip>
#include <iostream>
#include <stdexcept>

extern "C" int chps_quicked(const char*, int, const char*, int);

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
        const int score = chps_quicked(a.data(), static_cast<int>(a.size()), b.data(), static_cast<int>(b.size()));
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        if (score < 0) throw std::runtime_error("QuickEd failed to return a distance");
        rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
        usage.ru_maxrss /= 1024;
#endif
        // Some periodic sequences give a suboptimal upstream CIGAR. Keep the
        // comparison outside both the QuickEd timer and its RSS measurement.
        const auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
        auto oracle = edlibAlign(a.data(), static_cast<int>(a.size()), b.data(), static_cast<int>(b.size()), config);
        const int exact = oracle.editDistance;
        const auto status = oracle.status;
        edlibFreeAlignResult(oracle);
        if (status != EDLIB_STATUS_OK || exact < 0) throw std::runtime_error("edlib verification failed");
        if (score != exact)
            throw std::runtime_error("QuickEd returned " + std::to_string(score) + ", edlib returned " +
                                     std::to_string(exact));
        std::cout << std::setprecision(12) << "{\"algorithm\":\"quicked\",\"verified_by\":\"edlib\",\"m\":" << a.size()
                  << ",\"n\":" << b.size() << ",\"distance\":" << score << ",\"seconds\":" << elapsed.count()
                  << ",\"peak_rss_kib\":" << usage.ru_maxrss << "}\n";
    } catch (const std::exception& e) {
        std::cerr << "editdist-quicked: " << e.what() << '\n';
        return 1;
    }
}
