// Commands on stdin trigger one complete calculation on a resident pair.
#include "chps/banded.hpp"
#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include <sys/resource.h>
#include <algorithm>
#include <chrono>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>

using chps::Result;
using View = std::string_view;

Result align(View alg, View a, View b) {
#ifdef CHPS_REPEAT_INITIAL_ONLY
    if (alg == "nvpd-initial") return chps::nvpd_omp(a, b);
#else
    if (alg == "nvpd-simd") return chps::nvpd_omp(a, b);
    if (alg == "seq") return chps::sequential(a, b);
    if (alg == "diagonal") return chps::antidiagonal(a, b);
    if (alg == "edlib") return chps::edlib_distance(a, b);
    if (alg == "wfa2") return chps::wfa_distance(a, b);
    if (alg == "banded") return chps::banded(a, b);
    if (alg == "banded-trim") return chps::banded_trim(a, b);
    if (alg == "seq-trace") return chps::sequential(a, b, true);
    if (alg == "hirschberg-trace") return chps::hirschberg(a, b, true);
#ifdef CHPS_CUDA
    if (alg == "nvpd-cuda") return chps::nvpd_cuda(a, b);
#endif
#endif
    throw std::invalid_argument("unknown algorithm");
}

void replay(View a, View b, const Result& result) {
    size_t i = 0, j = 0;
    int count = 0, edits = 0;
    for (char op : result.cigar) {
        if (op >= '0' && op <= '9') {
            count = 10 * count + op - '0';
            continue;
        }
        if (count <= 0 || View("=XID").find(op) == View::npos) throw std::runtime_error("bad CIGAR");
        for (int k = 0; k < count; ++k) {
            if ((op != 'I' && i == a.size()) || (op != 'D' && j == b.size()))
                throw std::runtime_error("CIGAR out of range");
            if ((op == '=' || op == 'X') && (a[i] == b[j]) != (op == '='))
                throw std::runtime_error("CIGAR mismatch");
            i += op != 'I';
            j += op != 'D';
            edits += op != '=';
        }
        count = 0;
    }
    if (count || i != a.size() || j != b.size() || edits != result.distance)
        throw std::runtime_error("incorrect CIGAR");
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const auto records = chps::read_fasta(argv[1]);
    if (records.size() != 2) return 2;
    const View a = records[0].bases, b = records[1].bases, alg = argv[2];
    const int expected = std::stoi(argv[3]);
    std::cout << "ready" << std::endl;
    std::string command;
    while (std::getline(std::cin, command)) {
        if (command != "run") return 2;
        const auto begin = std::chrono::steady_clock::now();
        const Result result = align(alg, a, b);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        if (result.distance != expected) throw std::runtime_error("score differs from oracle");
        if (alg.ends_with("-trace")) replay(a, b, result);
        rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage)) return 2;
        auto rss = usage.ru_maxrss;
#ifdef __APPLE__
        rss /= 1024;
#endif
        std::cout << a.size() << ',' << b.size() << ',' << result.distance << ','
                  << std::setprecision(12) << seconds << ',' << rss << ','
                  << result.device_bytes << std::endl;
    }
}
