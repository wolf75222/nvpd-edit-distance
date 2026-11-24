// Linux benchmark of sequence decoding and requested C++ allocations.
#include "chps/fasta.hpp"
#include "../counters.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <new>
#include <string_view>
#include <sys/resource.h>

namespace {
bool counting = false;
std::uint64_t calls = 0, requested = 0;
}

void* operator new(std::size_t n) {
    if (void* ptr = std::malloc(n ? n : 1)) {
        if (counting) { ++calls; requested += n; }
        return ptr;
    }
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

int main(int argc, char** argv) {
#ifdef CHPS_PANEL_FAI
    if (argc != 3 && argc != 4) return 2;
    const std::string fasta = argv[1];
    const std::string index = argc == 4 ? argv[3] : "";
#else
    if (argc != 3) return 2;
#endif
    using Clock = std::chrono::steady_clock;
    const std::string_view mode = argv[2];
    if (mode != "reader" && mode != "profile" && mode != "check") return 2;
    std::vector<chps::Record> rows;
    for (int run = 0; run < (mode == "check" ? 1 : 2); ++run) {
        rows.clear();
        ReadCounters counters;
        const bool profile = mode == "profile" && run;
        if (profile) counters.start();
        const auto start = Clock::now();
        counting = profile;
#ifdef CHPS_PANEL_FAI
        rows = index.empty() ? chps::read_fasta(fasta) : chps::read_fasta(fasta, index);
#else
        rows = chps::read_fasta(argv[1]);
#endif
        counting = false;
        const double seconds = std::chrono::duration<double>(Clock::now()-start).count();
        if (profile) counters.stop();
        struct rusage rss{};
        if (getrusage(RUSAGE_SELF, &rss) || rows.size()!=2) return 1;
        std::cout << rows[0].bases.size() << ',' << rows[1].bases.size() << ",-1,"
                  << std::setprecision(12) << seconds << ',' << rss.ru_maxrss << '\n';
        if (profile) {
            counters.print();
            std::cout << "allocations," << calls << "\nrequested_bytes," << requested << '\n';
        }
        if (mode == "check")
            for (const auto& row : rows) std::cout << '>' << row.id << '\n' << row.bases << '\n';
    }
}
