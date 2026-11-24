#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/split.hpp"
#include <chrono>
#ifdef __linux__
#include "../../human/chromosome-panel/counters.hpp"
#endif
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) return 2;
    const auto left = chps::read_fasta(argv[1]), right = chps::read_fasta(argv[2]);
    if (left.size() != 1 || right.size() != 1) throw std::runtime_error("expected two single-record files");
    const auto& a = left.front().bases;
    const auto& b = right.front().bases;
    const std::string_view method = argv[3];
    for (int rep = 0; rep < 2; ++rep) {
#ifdef __linux__
        ReadCounters counters;
        const bool profile = argc == 5 && std::string_view(argv[4]) == "counters" && rep == 1;
        if (profile) counters.start();
#endif
        const auto start = std::chrono::steady_clock::now();
        const auto result = [&] {
            if (method == "edlib") return chps::edlib_distance(a,b);
            if (method == "wfa2") return chps::wfa_distance(a,b);
            if (method == "biwfa") return chps::biwfa(a,b);
            if (method == "split") return chps::nvpd_split(a,b);
            if (method == "sparse") return chps::nvpd_sparse(a,b);
            throw std::invalid_argument("unknown method");
        }();
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
#ifdef __linux__
        if (profile) counters.stop();
#endif
        rusage ru{};
        if (getrusage(RUSAGE_SELF,&ru)) return 3;
        std::cout << a.size() << ',' << b.size() << ',' << result.distance << ',' << std::setprecision(12)
                  << dt << ',' << ru.ru_maxrss << '\n';
#ifdef __linux__
        if (profile) counters.print();
#endif
    }
}
