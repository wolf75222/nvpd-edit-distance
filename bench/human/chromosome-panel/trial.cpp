// One process, one complete human chromosome pair. All times exclude output.
#include "counters.hpp"
#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <sys/resource.h>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    using Clock = std::chrono::steady_clock;
    const std::string_view method = argv[2];
    const bool reading = method == "reader" || method == "profile" || method == "check";
    if (!reading && method != "edlib" && method != "wfa2") return 2;
    std::vector<chps::Record> pair;
    if (!reading) pair = chps::read_fasta(argv[1]);
    long long previous = -1;
    const int runs = method == "check" ? 1 : 2;
    for (int run = 0; run < runs; ++run) {
        if (reading) pair.clear();
        ReadCounters counters;
        const bool profile = method == "profile" && run;
        if (profile) counters.start();
        const auto start = Clock::now();
        if (reading) pair = chps::read_fasta(argv[1]);
        if (pair.size() != 2) throw std::runtime_error("expected two FASTA records");
        const auto& a = pair[0].bases;
        const auto& b = pair[1].bases;
        const long long score = reading ? -1 : (method == "edlib" ? chps::edlib_distance(a,b) : chps::wfa_distance(a,b)).distance;
        const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
        if (profile) counters.stop();
        if (run && score != previous) throw std::runtime_error("warm and timed scores differ");
        previous = score;
        struct rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage)) throw std::runtime_error("getrusage failed");
        std::cout << a.size() << ',' << b.size() << ',' << score << ',' << std::setprecision(12)
                  << seconds << ',' << usage.ru_maxrss << '\n';
        if (profile) counters.print();
        if (method == "check") {
            for (const auto& row : pair)
                std::cout << '>' << row.id << '\n' << row.bases << '\n';
        }
    }
}
