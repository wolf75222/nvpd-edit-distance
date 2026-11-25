#include "chps/fasta.hpp"
#include "chps/split.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
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
#include <utility>
#include <vector>

#include <omp.h>
#ifdef __linux__
#include <sched.h>
#endif
#include <sys/mman.h>
#include <unistd.h>

namespace {
using Clock = std::chrono::steady_clock;
using Time = Clock::time_point;
constexpr int chromosomes = 4;
constexpr std::array<int, 2> cpus{0, 72};

double seconds(Time start, Time end) {
    return std::chrono::duration<double>(end - start).count();
}

void verify_cpu(int wanted) {
#ifdef __linux__
    cpu_set_t mask;
    if (sched_getaffinity(0, sizeof(mask), &mask) || CPU_COUNT(&mask) != 1 ||
        !CPU_ISSET(wanted, &mask) || sched_getcpu() != wanted)
        throw std::runtime_error("Worker not pinned to CPU " + std::to_string(wanted));
#else
    (void)wanted;
    throw std::runtime_error("CPU affinity verification requires Linux");
#endif
}

// Guard pages distinguish each replica in /proc/PID/numa_maps during the audit.
class Replica {
    void* allocation_;
    std::size_t mapped_;
    std::string_view sequence_;
public:
    explicit Replica(std::string_view original) {
        const auto page_size = sysconf(_SC_PAGESIZE);
        if (page_size <= 0) throw std::runtime_error("Unknown system page size");
        const auto page = static_cast<std::size_t>(page_size);
        if (original.empty() || original.size() > std::numeric_limits<std::size_t>::max() - 3 * page)
            throw std::length_error("Invalid chromosome length");
        const auto rounded = (original.size() + page - 1) / page * page;
        mapped_ = rounded + 2 * page;
        allocation_ = mmap(nullptr, mapped_, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (allocation_ == MAP_FAILED) throw std::system_error(errno, std::generic_category(), "mmap");
        auto* data = static_cast<char*>(allocation_) + page;
        if (mprotect(data, rounded, PROT_READ | PROT_WRITE)) {
            const auto error = errno;
            munmap(allocation_, mapped_);
            throw std::system_error(error, std::generic_category(), "mprotect");
        }
        std::memcpy(data, original.data(), original.size());
        sequence_ = {data, original.size()};
    }
    ~Replica() { munmap(allocation_, mapped_); }
    Replica(const Replica&) = delete;
    Replica& operator=(const Replica&) = delete;
    std::string_view view() const { return sequence_; }
};

struct Pair {
    std::string path;
    std::array<std::size_t, 2> lengths{};
    int score{};
    std::array<std::string, 2> originals;
    std::array<std::optional<Replica>, 2> replicas;

    std::string_view sequence(int side, bool local_copy) const {
        return local_copy ? replicas[side]->view() : std::string_view(originals[side]);
    }
};

using Panel = std::array<Pair, chromosomes>;

void manifest(const char* path, Panel& panel) {
    std::ifstream file(path);
    for (auto& pair : panel) {
        if (!(file >> pair.path >> pair.lengths[0] >> pair.lengths[1] >> pair.score))
            throw std::runtime_error("Expected four rows in chromosome manifest");
    }
    std::string extra;
    if (file >> extra) throw std::runtime_error("Unexpected extra row in manifest");
}

void read_panel(Panel& panel) {
    for (auto& pair : panel) {
        auto records = chps::read_fasta(pair.path);
        if (records.size() != 2) throw std::runtime_error("Expected two records in " + pair.path);
        for (int side = 0; side < 2; ++side) {
            if (records[side].bases.size() != pair.lengths[side])
                throw std::runtime_error("FASTA length differs from manifest");
            pair.originals[side] = std::move(records[side].bases);
        }
    }
}

void require_team(int thread) {
    if (omp_get_num_threads() != 2) throw std::runtime_error("Expected two OpenMP workers");
    verify_cpu(cpus[thread]);
}

void make_replicas(Panel& panel) {
    std::array<std::exception_ptr, 2> errors{};
    #pragma omp parallel num_threads(2)
    {
        const auto tid = omp_get_thread_num();
        try {
            require_team(tid);
            if (tid == 1) {
                for (auto& pair : panel)
                    for (int side = 0; side < 2; ++side)
                        pair.replicas[side].emplace(pair.originals[side]);
            }
        } catch (...) {
            errors[tid] = std::current_exception();
        }
    }
    for (auto error : errors) if (error) std::rethrow_exception(error);
}

struct Worker {
    std::vector<Time> stops;
    std::exception_ptr error;
    int scores{};
};

void align(const Panel& panel, bool local, bool prefix, int count, int rotation,
           std::array<Worker, 2>& workers) {
    const int rounds = (count + 1) / 2;
    #pragma omp parallel num_threads(2)
    {
        const int tid = omp_get_thread_num();
        auto& worker = workers[tid];
        try { require_team(tid); }
        catch (...) { worker.error = std::current_exception(); }
        for (int round = 0; round < rounds; ++round) {
            if (!worker.error && 2 * round + tid < count) {
                try {
                    const auto& pair = panel[(rotation + tid + round / chromosomes + round % chromosomes) % chromosomes];
                    const bool replica = local && tid == 1;
                    if (chps::nvpd_sparse(pair.sequence(0, replica), pair.sequence(1, replica)).distance != pair.score)
                        throw std::runtime_error("Score differs from oracle for " + pair.path);
                    ++worker.scores;
                    if (prefix) worker.stops[round] = Clock::now();
                    else if (2 * (round + 1) + tid >= count) worker.stops[0] = Clock::now();
                } catch (...) { worker.error = std::current_exception(); }
            }
            if (prefix) {
                #pragma omp barrier
            }
        }
    }
    for (const auto& worker : workers) if (worker.error) std::rethrow_exception(worker.error);
    if (workers[0].scores + workers[1].scores != count)
        throw std::runtime_error("Incomplete set of chromosome scores");
}

void audit(const Panel& panel, bool replicated) {
    for (int copy = 0; copy <= static_cast<int>(replicated); ++copy) {
        for (int chrom = 0; chrom < chromosomes; ++chrom) {
            for (int side = 0; side < 2; ++side) {
                const auto s = panel[chrom].sequence(side, copy == 1);
                std::cout << "buffer " << copy << ' ' << chrom << ' ' << side << ' ' << cpus[copy] << ' '
                          << reinterpret_cast<std::uintptr_t>(s.data()) << ' ' << s.size() << '\n';
            }
        }
    }
    std::cout << "audit" << std::endl;
    std::string ack;
    if (!std::getline(std::cin, ack) || ack != "done") throw std::runtime_error("Incomplete residency audit");
}
} // namespace

int main(int argc, char** argv) try {
    if ((argc != 6 && argc != 7) || (argc == 7 && std::string_view(argv[6]) != "--audit"))
        throw std::invalid_argument("Usage: grace_amortization manifest shared|local direct|prefix K rotation [--audit]");
    const std::string_view mode = argv[2], kind = argv[3];
    const bool local = mode == "local", prefix = kind == "prefix";
    const int count = std::stoi(argv[4]), rotation = std::stoi(argv[5]);
    if ((!local && mode != "shared") || (!prefix && kind != "direct") || count < 1 || count > 1000000 ||
        (prefix && count % 2) || rotation < 0 || rotation >= chromosomes)
        throw std::invalid_argument("Invalid configuration, score count or chromosome rotation");
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    omp_set_num_threads(1);
    verify_cpu(cpus[0]);
    Panel panel;
    manifest(argv[1], panel);
    std::array<Worker, 2> workers;
    for (auto& worker : workers) worker.stops.resize(prefix ? count / 2 : 1);

    const auto start = Clock::now();
    read_panel(panel);
    const auto loaded = Clock::now();
    if (local) make_replicas(panel);
    const auto ready = local ? Clock::now() : loaded;
    align(panel, local, prefix, count, rotation, workers);
    std::cout << std::setprecision(12);
    for (std::size_t i = 0; i < workers[0].stops.size(); ++i) {
        const auto stop = std::max(workers[0].stops[i], workers[1].stops[i]);
        const auto scored = prefix ? 2 * static_cast<int>(i + 1) : count;
        std::cout << "timing " << scored << ' ' << seconds(start, stop) << ' ' << seconds(start, loaded)
                  << ' ' << seconds(loaded, ready) << ' ' << seconds(start, ready) << ' '
                  << seconds(ready, stop) << ' ' << scored << '\n';
    }
    if (argc == 7) audit(panel, local);
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
