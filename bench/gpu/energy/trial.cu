#include "chps/cuda.hpp"
#include "chps/fasta.hpp"
#include <cuda_runtime.h>
#include <nvml.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

using Clock = std::chrono::steady_clock;
using Seconds = std::chrono::duration<double>;

void nvcheck(nvmlReturn_t status) {
    if (status != NVML_SUCCESS) throw std::runtime_error(nvmlErrorString(status));
}

// The cumulative counter is preferred. Devices without it use trapezoidal
// integration of the driver power sensor at its native update cadence.
class DeviceMeter {
    nvmlDevice_t device_{};
    unsigned long long start_mj_{};
    std::atomic<bool> done_{false};
    std::thread thread_;
    double integrated_j_ = 0;
    nvmlReturn_t poll_error_ = NVML_SUCCESS;
    bool counter_ = false;

public:
    DeviceMeter() {
        char pci[32]{};
        if (cudaDeviceGetPCIBusId(pci, sizeof(pci), 0) != cudaSuccess)
            throw std::runtime_error("no PCI address for CUDA device 0");
        nvcheck(nvmlInit_v2());
        nvcheck(nvmlDeviceGetHandleByPciBusId_v2(pci, &device_));
        counter_ = nvmlDeviceGetTotalEnergyConsumption(device_, &start_mj_) == NVML_SUCCESS;
    }
    ~DeviceMeter() { done_ = true; if (thread_.joinable()) thread_.join(); nvmlShutdown(); }
    const char* sensor() const { return counter_ ? "total_energy_mJ" : "power_trapezoid"; }
    void start() {
        if (counter_) { nvcheck(nvmlDeviceGetTotalEnergyConsumption(device_, &start_mj_)); return; }
        unsigned initial;
        nvcheck(nvmlDeviceGetPowerUsage(device_, &initial));
        thread_ = std::thread([this, initial] {
            auto previous = Clock::now();
            unsigned last = initial;
            do {
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
                unsigned power;
                auto status = nvmlDeviceGetPowerUsage(device_, &power);
                if (status != NVML_SUCCESS) { poll_error_ = status; return; }
                auto now = Clock::now();
                integrated_j_ += Seconds(now - previous).count() * (last + power) / 2000.;
                previous = now;
                last = power;
            } while (!done_);
        });
    }
    double stop() {
        if (counter_) {
            unsigned long long end;
            nvcheck(nvmlDeviceGetTotalEnergyConsumption(device_, &end));
            if (end <= start_mj_) throw std::runtime_error("nonincreasing energy counter");
            return (end - start_mj_) / 1000.;
        }
        done_ = true;
        thread_.join();
        nvcheck(poll_error_);
        return integrated_j_;
    }
};

int main(int argc, char** argv) try {
    if (argc != 4) throw std::runtime_error("trial FASTA expected_distance minimum_seconds");
    const auto seq = chps::read_fasta(argv[1]);
    if (seq.size() != 2) throw std::runtime_error("two records required");
    const int expected = std::stoi(argv[2]);
    const double duration = std::stod(argv[3]);
    chps::CudaWorkspace workspace(chps::CudaHalo::local);
    auto solve = [&] {
        const int score = workspace.distance(seq[0].bases, seq[1].bases).distance;
        if (score != expected) throw std::runtime_error("mismatched distance");
    };
    const auto cold = Clock::now();
    solve();
    const double cold_s = Seconds(Clock::now() - cold).count();
    if (duration == 0.) return 0;  // Single solve for Nsight collection.
    const auto warm = Clock::now();
    solve();
    const double warm_s = Seconds(Clock::now() - warm).count();
    const int iterations = std::max(1, static_cast<int>(std::ceil(duration / warm_s)));
    DeviceMeter meter;
    meter.start();
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) solve();
    const double seconds = Seconds(Clock::now() - start).count();
    const double joules = meter.stop();
    std::cout << std::setprecision(12) << expected << ',' << seq[0].bases.size() << ','
              << seq[1].bases.size() << ',' << cold_s << ',' << warm_s << ',' << iterations << ','
              << seconds << ',' << joules << ',' << meter.sensor() << '\n';
    return 0;
} catch (const std::exception& ex) {
    std::cerr << ex.what() << '\n';
    return 1;
}
