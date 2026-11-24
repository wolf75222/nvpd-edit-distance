#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/split.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <system_error>
#include <vector>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

using Clock = std::chrono::steady_clock;
using View = std::string_view;
constexpr std::size_t huge = 2 * 1024 * 1024;
constexpr int calls = 3;

void checked(bool ok, const char* action) {
    if (!ok) throw std::system_error(errno, std::generic_category(), action);
}

class InputMapping {
    void* reservation_;
    std::size_t reserved_;
    char* body_;
    std::size_t bytes_;
    std::size_t offset_;
    static std::size_t round_up(std::size_t size) { return (size + huge - 1) / huge * huge; }
public:
    InputMapping(View a, View b, bool thp) {
        offset_ = round_up(a.size());
        bytes_ = offset_ + round_up(b.size());
        reserved_ = bytes_ + 2 * huge;
        reservation_ = mmap(nullptr, reserved_, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        checked(reservation_ != MAP_FAILED, "reserve inputs");
        auto start = reinterpret_cast<std::uintptr_t>(reservation_);
        body_ = reinterpret_cast<char*>(round_up(start + 1));
        checked(mprotect(body_, bytes_, PROT_READ | PROT_WRITE) == 0, "protect inputs");
        checked(madvise(body_, bytes_, thp ? MADV_HUGEPAGE : MADV_NOHUGEPAGE) == 0, "madvise inputs");
        std::memcpy(body_, a.data(), a.size());
        std::memset(body_ + a.size(), 0, offset_ - a.size());
        std::memcpy(body_ + offset_, b.data(), b.size());
        std::memset(body_ + offset_ + b.size(), 0, bytes_ - offset_ - b.size());
    }
    InputMapping(const InputMapping&) = delete;
    InputMapping& operator=(const InputMapping&) = delete;
    ~InputMapping() { munmap(reservation_, reserved_); }
    View first(std::size_t n) const { return {body_, n}; }
    View second(std::size_t n) const { return {body_ + offset_, n}; }
    std::size_t capacity() const { return bytes_; }

    std::size_t huge_bytes() const {
        std::ifstream smaps("/proc/self/smaps");
        bool in_region = false;
        for (std::string line; std::getline(smaps, line);) {
            unsigned long low, high, kib;
            if (std::sscanf(line.c_str(), "%lx-%lx", &low, &high) == 2)
                in_region = low == reinterpret_cast<std::uintptr_t>(body_) && high - low == bytes_;
            if (in_region && std::sscanf(line.c_str(), "AnonHugePages: %lu kB", &kib) == 1)
                return kib * 1024;
        }
        throw std::runtime_error("input VMA not found in smaps");
    }
};

struct Event { unsigned type; std::uint64_t code; };
struct CounterGroup {
    std::vector<int> descriptors;
    explicit CounterGroup(std::initializer_list<Event> events) {
        for (auto event : events) {
            perf_event_attr attr{};
            attr.size = sizeof(attr);
            attr.type = event.type;
            attr.config = event.code;
            attr.disabled = attr.exclude_kernel = attr.exclude_hv = 1;
            attr.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
            int leader = descriptors.empty() ? -1 : descriptors.front();
            int fd = syscall(SYS_perf_event_open, &attr, 0, -1, leader, 0);
            checked(fd >= 0, "open user PMU group");
            descriptors.push_back(fd);
        }
    }
    CounterGroup(const CounterGroup&) = delete;
    CounterGroup& operator=(const CounterGroup&) = delete;
    ~CounterGroup() { for (auto fd : descriptors) close(fd); }
    void enable() const { checked(ioctl(descriptors[0], PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == 0, "enable counters"); }
    void finish() const {
        checked(ioctl(descriptors[0], PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP) == 0, "stop counters");
        std::vector<std::uint64_t> values(3 + descriptors.size());
        auto bytes = static_cast<ssize_t>(values.size() * sizeof(values[0]));
        checked(read(descriptors[0], values.data(), bytes) == bytes, "read counters");
        if (values[0] != descriptors.size() || !values[1] || values[2] < .99 * values[1])
            throw std::runtime_error("PMU group has less than 99% active time");
        std::cout << ',' << double(values[2])/values[1];
        for (std::size_t i = 3; i < values.size(); ++i) std::cout << ',' << values[i];
    }
};

double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void measure_reuse(View a, View b, int oracle) {
    std::array<std::pair<int, bool>, 6> configs{{{1,false},{1,true},{192,false},{192,true},{384,false},{384,true}}};
    std::mt19937 rng(601);
    for (int block = 1; block <= 5; ++block) {
        std::shuffle(configs.begin(), configs.end(), rng);
        for (const auto& [repeats, thp] : configs) {
            if (chps::nvpd_split(a, b).distance != oracle) throw std::runtime_error("warm-up score differs");
            auto start = Clock::now();
            InputMapping input(a, b, thp);
            auto x = input.first(a.size()), y = input.second(b.size());
            auto mapped = Clock::now();
            for (int k = 0; k < repeats; ++k)
                if (chps::nvpd_split(x, y).distance != oracle) throw std::runtime_error("reused score differs");
            const double total = elapsed(start);
            std::cout << block << ',' << (thp ? "thp" : "base") << ',' << repeats << ',' << oracle << ','
                      << std::chrono::duration<double>(mapped-start).count() << ',' << total << ','
                      << input.capacity() << ',' << input.huge_bytes() << '\n' << std::flush;
        }
    }
}

int main(int argc, char** argv) try {
    bool reuse = argc == 5 && View(argv[4]) == "reuse";
    if (argc != 4 && !reuse) throw std::runtime_error("usage: trial CHROMOSOME FASTA EXPECTED_DISTANCE [reuse]");
    auto records = chps::read_fasta(argv[2]);
    if (records.size() != 2) throw std::runtime_error("a pair of FASTA records is required");
    auto &a = records[0].bases, &b = records[1].bases;
    const int expected = std::stoi(argv[3]);
    const auto score = [](View x, View y, bool wfa) {
        return (wfa ? chps::wfa_distance(x, y) : chps::nvpd_split(x, y)).distance;
    };
    const int oracle = chps::edlib_distance(a, b).distance;
    if (oracle != expected) throw std::runtime_error("edlib disagrees with expected score");
    if (reuse) {
        std::cout << std::setprecision(11);
        measure_reuse(a, b, oracle);
        return 0;
    }
    std::array<int, 4> configs{0, 1, 2, 3};
    std::mt19937 rng(8703);
    std::cout << std::setprecision(11);
    for (int block = 1; block <= 9; ++block) {
        std::shuffle(configs.begin(), configs.end(), rng);
        for (int config : configs) {
            bool thp = (config & 1), wfa = (config & 2);
            auto start = Clock::now();
            InputMapping input(a, b, thp);
            const double setup = elapsed(start);
            View x = input.first(a.size()), y = input.second(b.size());
            if (score(x, y, wfa) != oracle) throw std::runtime_error("warm-up score differs");
            std::cout << argv[1] << ',' << block << ',' << (wfa ? "wfa2" : "split") << ','
                      << (thp ? "thp" : "base") << ',' << a.size() << ',' << b.size() << ',' << expected
                      << ',' << input.capacity() << ',' << input.huge_bytes() << ',' << setup << ',' << calls;
            const auto batch = [&](const CounterGroup& group) {
                group.enable();
                auto begin = Clock::now();
                for (int i = 0; i < calls; ++i)
                    if (score(x, y, wfa) != oracle) throw std::runtime_error("measured score differs");
                auto seconds = elapsed(begin);
                group.finish();
                std::cout << ',' << seconds;
            };
            batch(CounterGroup{{PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
                               {PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
                               {PERF_TYPE_HW_CACHE, PERF_COUNT_HW_CACHE_L1D | (PERF_COUNT_HW_CACHE_RESULT_MISS << 16)},
                               {PERF_TYPE_RAW, 0xf045}});
            // Zen 4 core refills from local or remote DRAM/MMIO, including prefetch.
            batch(CounterGroup{{PERF_TYPE_RAW, 0x4844}, {PERF_TYPE_RAW, 0x4843}, {PERF_TYPE_RAW, 0x485a}});
            std::cout << '\n' << std::flush;
        }
    }
} catch (const std::exception& err) {
    std::cerr << err.what() << '\n';
    return 1;
}
