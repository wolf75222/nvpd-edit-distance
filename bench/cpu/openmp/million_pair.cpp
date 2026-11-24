// Profile one cold, cost-only call, sampling worker placement after the timer.
#include "chps/dp.hpp"
#include "chps/fasta.hpp"
#include <omp.h>
#include <sched.h>
#include <chrono>
#include <iostream>
#include <set>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const auto records = chps::read_fasta(argv[1]);
    if (records.size() != 2) return 2;
    const auto& a = records[0].bases;
    const auto& b = records[1].bases;
    const auto start = std::chrono::steady_clock::now();
    const auto result = chps::nvpd_omp(a,b,false);
    const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::vector<int> cores(omp_get_max_threads(), -1);
    #pragma omp parallel
    cores[omp_get_thread_num()] = sched_getcpu();
    std::set<int> distinct(cores.begin(), cores.end());
    if (distinct.size()!=8 || distinct.count(-1)) return 1;
    std::cout << a.size() << ',' << b.size() << ',' << result.distance << ',' << seconds << ',';
    for (int core : distinct) std::cout << (core==*distinct.begin() ? "" : ";") << core;
    std::cout << '\n';
}
