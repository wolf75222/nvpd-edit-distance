#include "chps/adaptive.hpp"
#include "chps/banded.hpp"
#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "nvpd_band_row.hpp"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <omp.h>

class Gate {
    std::ofstream ctl;
    std::ifstream ack;
public:
    Gate() {
        if (const char* fifo = std::getenv("CHPS_PERF_CTL")) {
            ctl.open(fifo);
            ack.open(std::getenv("CHPS_PERF_ACK"));
            if (!ctl || !ack) throw std::runtime_error("PMU FIFO unavailable");
        }
    }
    void set(const char* state) {
        if (!ctl.is_open()) return;
        ctl << state << std::endl;
        std::string line;
        while (ack.peek() == 0) ack.get();
        if (!std::getline(ack,line) || line != "ack") throw std::runtime_error("PMU acknowledgement");
    }
};

template<class F> double timed(Gate& gate, F run) {
    gate.set("enable");
    const auto start = std::chrono::steady_clock::now();
    run();
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    gate.set("disable");
    return elapsed.count();
}

void alignment(const std::string& path, std::string_view mode) {
    const auto pair=chps::read_fasta(path);
    if (pair.size()!=2) throw std::runtime_error("expected two FASTA records");
    const auto& a=pair[0].bases;
    const auto& b=pair[1].bases;
    chps::AdaptiveWork work;
    auto run = [&]() {
        if (mode=="nvpd-adaptive") return chps::nvpd_adaptive_score(a,b,work);
        if (mode=="nvpd-omp") return chps::nvpd_omp(a,b).distance;
        if (mode=="banded-trim") return chps::banded_trim(a,b).distance;
        if (mode=="biwfa") return chps::biwfa(a,b).distance;
        if (mode=="wfa2") return chps::wfa_distance(a,b).distance;
        if (mode=="edlib") return chps::edlib_distance(a,b).distance;
        throw std::invalid_argument("unknown algorithm");
    };
    const int warm=run();
    Gate gate;
    int distance;
    const double seconds=timed(gate,[&]{distance=run();});
    if (distance!=warm) throw std::runtime_error("warm score mismatch");
    if (mode=="nvpd-omp") work.cells=static_cast<std::uint64_t>(a.size())*b.size();
    std::cout<<seconds<<','<<distance<<','<<work.cells<<','<<work.passes<<','<<work.threshold
             <<','<<work.prefix<<','<<work.suffix<<'\n';
}

// Memory clobber keeps all calibrated traversals visible at -O3 without
// serializing memory accesses in their inner loops (GCC/Clang benchmark only).
template<class T> inline void consume(const T* address) {
    asm volatile("" : : "g"(address) : "memory");
}

// Use a direct pointer so GCC emits contiguous vector loads. The modular
// unsigned sum keeps the printed checksum observable and defined.
__attribute__((noinline)) unsigned stream_read(const unsigned* data, int count) {
    unsigned sum=0;
    #pragma omp simd reduction(+:sum)
    for (int j=0;j<count;++j) sum+=data[j];
    return sum;
}

void calibration(std::string_view kind, int footprint_kib) {
    const int arrays=kind=="read"?1:3;
    const int count=footprint_kib*1024/(4*arrays);
    const int rounds=std::max(3,static_cast<int>((256ULL<<20)/(count*4ULL*arrays)));
    std::vector<unsigned> input(count), multiplier, output;
    std::iota(input.begin(),input.end(),0U);
    if (kind=="triad") { multiplier=input; output.resize(count); }
    std::vector<int> prev, curr, last;
    if (kind=="row") {
        prev.resize(count+3); curr.resize(count+3); last.resize(count+3);
        for (int j=1;j<=count;++j) { prev[j]=j%31; last[j]=j-j%4; }
    }
    int row_best=0;
    std::uint64_t checksum=0;
    auto sweep = [&] {
        if (kind=="read") {
            consume(input.data());
            checksum+=stream_read(input.data(),count);
        } else if (kind=="triad") {
            #pragma omp simd
            for (int j=0;j<count;++j) output[j]=input[j]+3U*multiplier[j];
            consume(output.data());
            checksum+=output[count/2];
        } else if (kind=="row") {
            row_best+=chps::detail::nvpd_band_row(prev.data(),curr.data(),last.data(),0,count,0,65535);
            consume(curr.data());
            checksum+=curr[count/2];
        } else throw std::invalid_argument("unknown calibration");
    };
    for (int r=0;r<rounds;++r) sweep();
    Gate gate;
    const double elapsed=timed(gate,[&]{for(int r=0;r<rounds;++r) sweep();});
    std::cout<<elapsed<<','<<count<<','<<rounds<<','<<checksum+row_best<<'\n';
}
int main(int argc,char** argv) {
    if(argc!=4) return 2;
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    std::cout<<std::setprecision(12);
    if(std::string_view(argv[1])=="align") alignment(argv[2],argv[3]);
    else if(std::string_view(argv[1])=="calibrate") calibration(argv[2],std::stoi(argv[3]));
    else return 2;
}
