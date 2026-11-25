#include "chps/fasta.hpp"
#include "chps/split.hpp"

#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include <omp.h>
#ifdef __linux__
#include <sched.h>
#endif
#include <sys/mman.h>
#include <unistd.h>

namespace {
using Clock = std::chrono::steady_clock;
constexpr int chromosomes = 4;
constexpr int strong_passes = 16;

class Sequence {
    char* allocation_;
    std::size_t total_;
    std::string_view sequence_;
public:
    // Guard pages keep this individual sequence identifiable in /proc/PID/numa_maps.
    explicit Sequence(std::string_view source) {
        const long page_size = sysconf(_SC_PAGESIZE);
        if (page_size <= 0) throw std::runtime_error("Cannot determine page size");
        const auto page = static_cast<std::size_t>(page_size);
        if (source.empty() || source.size() > std::numeric_limits<std::size_t>::max() - 3 * page)
            throw std::length_error("Invalid sequence length");
        const auto rounded = ((source.size() + page - 1) / page) * page;
        total_ = rounded + 2 * page;
        void* base = mmap(nullptr, total_, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (base == MAP_FAILED) throw std::system_error(errno, std::generic_category(), "mmap");
        allocation_ = static_cast<char*>(base);
        char* data = allocation_ + page;
        if (mprotect(data, rounded, PROT_READ | PROT_WRITE) < 0) {
            const int error = errno;
            munmap(base, total_);
            throw std::system_error(error, std::generic_category(), "mprotect");
        }
        std::memcpy(data, source.data(), source.size()); // First write from the bound owner.
        sequence_ = {data, source.size()};
    }
    ~Sequence() { munmap(allocation_, total_); }
    Sequence(const Sequence&) = delete;
    Sequence& operator=(const Sequence&) = delete;
    std::string_view view() const { return sequence_; }
};

struct Pair {
    std::array<std::optional<Sequence>, 2> bases;
    int score{};
};
using Panel = std::array<Pair, chromosomes>;

void check_cpu(int expected) {
#ifdef __linux__
    cpu_set_t mask;
    if (expected < 0 || expected >= CPU_SETSIZE || sched_getaffinity(0, sizeof(mask), &mask) != 0 ||
        CPU_COUNT(&mask) != 1 || !CPU_ISSET(expected, &mask) || sched_getcpu() != expected)
        throw std::runtime_error("Thread is not bound to requested CPU " + std::to_string(expected));
#else
    (void)expected;
    throw std::runtime_error("Measuring verified per-thread CPU affinity requires Linux");
#endif
}

template<class F>
double with_workers(int workers, const std::array<int, 2>& cpus, F work) {
    std::array<std::exception_ptr, 2> errors{};
    const auto start = Clock::now();
    #pragma omp parallel num_threads(workers)
    {
        const int tid = omp_get_thread_num();
        try {
            omp_set_num_threads(1);
            if (omp_get_num_threads() != workers) throw std::runtime_error("Unexpected OpenMP team size");
            check_cpu(cpus[tid]);
            work(tid);
        } catch (...) {
            errors[tid] = std::current_exception();
        }
    }
    const auto stop = Clock::now();
    for (const auto& error : errors)
        if (error) std::rethrow_exception(error);
    return std::chrono::duration<double>(stop - start).count();
}

void wait_for(std::string_view expected) {
    std::string line;
    if (!std::getline(std::cin, line) || line != expected)
        throw std::runtime_error("Expected control message " + std::string(expected));
}
} // namespace

int main(int argc, char** argv) try {
    if (argc != 5) throw std::invalid_argument("Usage: grace_numa inputs.tsv MODE cpu0 cpu1");
    const std::string_view mode = argv[2];
    const bool replicated = mode == "local" || mode == "remote" || mode == "twins0" || mode == "twins1";
    if (!replicated && mode != "single" && mode != "shared0" && mode != "shared1")
        throw std::invalid_argument("Unknown memory configuration");
    const int workers = mode == "single" ? 1 : 2;
    const int copies = replicated ? 2 : 1;
    const std::array<int, 2> cpus{std::stoi(argv[3]), std::stoi(argv[4])};
    if (cpus[0] == cpus[1]) throw std::invalid_argument("Expected two distinct CPU numbers");
    const auto owner = [mode](int copy) {
        if (mode == "remote") return 1 - copy;
        if (mode == "twins0") return 0;
        if (mode == "shared1" || mode == "twins1") return 1;
        return copy;
    };
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    omp_set_num_threads(1);
    std::ifstream input(argv[1]);
    std::array<Panel, 2> panels;
    std::uint64_t bases = 0;
    for (int k = 0; k < chromosomes; ++k) {
        std::string path;
        std::size_t m, n;
        int score;
        if (!(input >> path >> m >> n >> score)) throw std::runtime_error("Expected four chromosome pairs");
        const auto records = chps::read_fasta(path);
        if (records.size() != 2 || records[0].bases.size() != m || records[1].bases.size() != n)
            throw std::runtime_error("FASTA differs from manifest: " + path);
        bases += m + n;
        with_workers(workers, cpus, [&](int tid) {
            for (int copy = 0; copy < copies; ++copy) {
                if (owner(copy) != tid) continue;
                auto& pair = panels[copy][k];
                pair.score = score;
                for (int side = 0; side < 2; ++side) pair.bases[side].emplace(records[side].bases);
            }
        });
    }
    std::string extra;
    if (input >> extra) throw std::runtime_error("Surplus chromosome pairs");
    const auto measure = [&](int passes) {
        return with_workers(workers, cpus, [&](int tid) {
            const auto& panel = panels[replicated ? tid : 0];
            for (int pass = 0; pass < passes; ++pass) {
                for (int k = 0; k < chromosomes; ++k) {
                    const auto& pair = panel[(tid + pass + k) % chromosomes];
                    if (chps::nvpd_sparse(pair.bases[0]->view(), pair.bases[1]->view()).distance != pair.score)
                        throw std::runtime_error("NvPD score differs from the chromosome oracle");
                }
            }
        });
    };
    const double warmup = measure(1);
    for (int copy = 0; copy < copies; ++copy) {
        for (int k = 0; k < chromosomes; ++k) {
            for (int side = 0; side < 2; ++side) {
                const auto seq = panels[copy][k].bases[side]->view();
                std::cout << "buffer " << copy << ' ' << k << ' ' << side << ' ' << cpus[owner(copy)] << ' '
                          << reinterpret_cast<std::uintptr_t>(seq.data()) << ' ' << seq.size() << '\n';
            }
        }
    }
    std::cout << std::setprecision(12) << "ready " << workers << ' ' << strong_passes * chromosomes << ' '
              << strong_passes * bases << ' ' << warmup << std::endl;
    wait_for("run");
    const double elapsed = measure(strong_passes / workers);
    std::cout << "result " << elapsed << std::endl;
    wait_for("done");
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
