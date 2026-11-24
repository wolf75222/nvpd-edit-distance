#pragma once
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>

class ReadCounters {
    struct Event { const char* name; unsigned type; std::uint64_t code; };
    static constexpr std::array<Event,4> events{{
        {"instructions", PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
        {"cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
        {"L1-dcache-loads", PERF_TYPE_HW_CACHE, PERF_COUNT_HW_CACHE_L1D},
        {"L1-dcache-load-misses", PERF_TYPE_HW_CACHE,
         PERF_COUNT_HW_CACHE_L1D | (PERF_COUNT_HW_CACHE_RESULT_MISS << 16)}
    }};
    std::array<int,4> descriptors{{-1,-1,-1,-1}};

public:
    void start() {
        for (std::size_t i = 0; i < events.size(); ++i) {
            perf_event_attr attr{};
            attr.size = sizeof(attr);
            attr.type = events[i].type;
            attr.config = events[i].code;
            attr.disabled = 1;
            attr.exclude_kernel = attr.exclude_hv = 1;
            attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
            descriptors[i] = syscall(SYS_perf_event_open, &attr, 0, -1, descriptors[0], 0);
            if (descriptors[i] < 0) throw std::system_error(errno, std::generic_category(), events[i].name);
        }
        if (ioctl(descriptors[0], PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) < 0)
            throw std::system_error(errno, std::generic_category(), "enable PMU group");
    }
    void stop() {
        if (ioctl(descriptors[0], PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP) < 0)
            throw std::system_error(errno, std::generic_category(), "disable PMU group");
    }
    void print() const {
        for (std::size_t i = 0; i < events.size(); ++i) {
            std::array<std::uint64_t,3> data{};
            if (read(descriptors[i], data.data(), sizeof(data)) != sizeof(data) || !data[2])
                throw std::runtime_error("PMU group was not scheduled");
            std::cout << events[i].name << ',' << data[0] << ',' << 100.0 * data[2] / data[1] << '\n';
        }
    }
    ~ReadCounters() { for (int fd : descriptors) if (fd >= 0) close(fd); }
};
