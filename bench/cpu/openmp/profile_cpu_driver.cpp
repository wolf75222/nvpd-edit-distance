// Isolate one NvPD call from input generation, validation and warm-up.
#include "chps/dp.hpp"
#include <omp.h>
#include <sched.h>
#include <sys/resource.h>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 6) return 2;
    const std::string algorithm(argv[1]);
    const int m = std::stoi(argv[2]), n = std::stoi(argv[3]);
    const int threads = std::stoi(argv[4]), warmups = std::stoi(argv[5]);
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
    std::mt19937 rng(2026);
    constexpr char dna[] = "ACGT";
    std::string a(m, 'A'), b(n, 'A');
    for (auto& c : a) c = dna[rng() % 4];
    for (int j = 0; j < n; ++j) {
        b[j] = j < m ? a[j] : dna[rng() % 4];
        if (j < m && rng() % 20 == 0) {
            char c = b[j];
            while (c == b[j]) c = dna[rng() % 4];
            b[j] = c;
        }
    }
    std::uint64_t hash = 14695981039346656037ull;
    for (const auto& input : {&a, &b})
        for (unsigned char c : *input) hash = (hash ^ c) * 1099511628211ull;
    int max_gap = 0;
    for (char c : std::string_view("ACGT")) {
        int gap = 0;
        if (a.find(c) == std::string::npos) continue;
        for (char b_j : b) {
            gap = b_j == c ? 0 : gap + 1;
            if (gap > max_gap) max_gap = gap;
        }
    }
    const auto fn = algorithm == "seq" ? chps::sequential : chps::nvpd_omp;
    for (int warm = 0; warm < warmups; ++warm) fn(a, b, false);
    const auto* ctl = std::getenv("CHPS_PERF_CTL");
    std::ofstream control;
    std::ifstream ack;
    if (ctl) {
        control.open(ctl);
        ack.open(std::getenv("CHPS_PERF_ACK"));
        if (!control || !ack) throw std::runtime_error("perf control FIFO unavailable");
    }
    auto gate = [&](const char* op) {
        if (!ctl) return;
        control << op << std::endl;
        std::string line;
        // perf 6.9 sends a NUL after each newline-terminated acknowledgement.
        while (ack.peek() == 0) ack.get();
        if (!std::getline(ack, line) || line != "ack")
            throw std::runtime_error("perf did not acknowledge control command");
    };
    gate("enable");
    const auto start = std::chrono::steady_clock::now();
    const auto result = fn(a, b, false);
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    gate("disable");
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::vector<int> cores(threads);
    #pragma omp parallel
    cores[omp_get_thread_num()] = sched_getcpu();
    std::cout << result.distance << ',' << std::setprecision(12) << elapsed.count()
              << ',' << usage.ru_maxrss << ',' << std::hex << hash << std::dec << ',';
    for (int t = 0; t < threads; ++t) std::cout << (t ? ";" : "") << cores[t];
    std::cout << ',' << max_gap << '\n';
}
