// Enumerate the source-level global word addresses of cost-only NvPD on each warp.
// Geometry and predicates mirror src/cuda.cu. No counters or timing estimates.
#include "chps/fasta.hpp"
#include "dp_internal.hpp"
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
constexpr int warp = 32, block = 256, tile_rows = 8;
enum Access { mit, diagonal, upper, gather, output, kinds };
constexpr const char* names[kinds] = {"mit", "diagonal", "upper", "gather", "store"};

struct Count {
    std::uint64_t requests{}, lanes{}, words{}, sectors{};
};

struct WarpAccess {
    Count count;
    int last_word = -1;
    void add(int word) {
        // Each individual instruction's addresses are monotonically nondecreasing:
        // this includes k(j)-1 because the last matching position never decreases.
        if (word < last_word) throw std::logic_error("non-monotone warp addresses");
        if (last_word == -1 || word/8 != last_word/8) ++count.sectors;
        if (word != last_word) ++count.words;
        ++count.lanes;
        last_word = word;
    }
    void finish(Count& total) const {
        total.requests += count.lanes != 0;
        total.lanes += count.lanes;
        total.words += count.words;
        total.sectors += count.sectors;
    }
};

void check_model() {
    WarpAccess aligned, shifted, broadcast;
    for (int j=0; j<warp; ++j) {
        aligned.add(j);
        shifted.add(j+1);
        broadcast.add(30);
    }
    Count sum{};
    broadcast.finish(sum);
    if (aligned.count.sectors != 4 || shifted.count.sectors != 5 || sum.requests != 1 ||
        sum.sectors != 1 || sum.words != 1 || sum.lanes != warp)
        throw std::logic_error("warp reference pattern failed");
}

void check_cuda(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(cudaGetErrorString(status));
}

void model(const std::string& label, const std::string& a, const std::string& b,
           const chps::detail::MatchIndex& index, int reach, int sms, int score, bool adaptive) {
    const int m = a.size(), n = b.size(), stride = n+1;
    const bool tiled = adaptive && m >= tile_rows && reach <= 64;
    const int height = tiled ? tile_rows : 1;
    const int width = tiled && stride / (4*block) >= 2*sms ? 4*block : block;
    int launches = 0;
    std::array<Count, kinds> totals{};
    for (int first=0; first<m; first+=height, ++launches) {
        const int rows = std::min(height, m-first);
        const int input = (launches%2)*stride, dest = (1-launches%2)*stride;
        for (int begin=0; begin<=n; begin+=width) {
            const int end = std::min(begin+width-1, n);
            for (int t=0; t<rows; ++t) {
                const int left = std::max(0, begin-(rows-t-1)*reach);
                const int symbol = index.symbol[static_cast<unsigned char>(a[first+t])];
                const int base = symbol*stride;
                // Consecutive 32-lane groups stay distinct across all 256-thread loop passes.
                for (int start=left; start<=end; start+=warp) {
                    std::array<WarpAccess, kinds> accesses{};
                    for (int j=start; j<=std::min(start+warp-1,end); ++j) {
                        if (t+1 == rows) accesses[output].add(dest+j);
                        if (j == 0) continue;
                        accesses[mit].add(base+j);
                        if (t > 0) continue; // The previous DP row is in shared memory.
                        const int k = index.entries[base+j];
                        accesses[diagonal].add(input+j-1);
                        if (k != j) {
                            accesses[upper].add(input+j);
                            if (k) accesses[gather].add(input+k-1);
                        }
                    }
                    for (int s=0; s<kinds; ++s) accesses[s].finish(totals[s]);
                }
            }
        }
    }
    if (totals[output].lanes != static_cast<std::uint64_t>(stride)*launches ||
        (!tiled && totals[mit].lanes != static_cast<std::uint64_t>(m)*n))
        throw std::logic_error("modeled cells or output stores differ from DP dimensions");
    for (int s=0; s<kinds; ++s) {
        const auto& c = totals[s];
        std::cout << label << ',' << (adaptive ? "adaptive" : "single_row") << ','
                  << (tiled ? "row_tile" : "next_row") << ',' << m << ',' << n << ',' << reach
                  << ',' << sms << ',' << width << ',' << launches << ',' << score << ',' << names[s]
                  << ',' << c.requests << ',' << c.lanes << ',' << c.words << ',' << c.sectors << '\n';
    }
}
}

int main(int argc, char** argv) try {
    if (argc != 4) throw std::invalid_argument("gpu_sectors {model|capture} label pair.fa");
    auto data = chps::read_fasta(argv[3]);
    if (data.size() != 2) throw std::invalid_argument("expected two records");
    auto& a = data[0].bases;
    auto& b = data[1].bases;
    if (a.size() > b.size()) a.swap(b);
    if (a.empty() || b.size() > 1'000'000) throw std::invalid_argument("model supports 1..1M nt");
    const int expected = chps::sequential(a,b).distance;
    if (chps::nvpd_cuda(a,b).distance != expected) throw std::runtime_error("CUDA score differs from scalar");
    if (std::string(argv[1]) == "capture") {
        check_cuda(cudaProfilerStart());
        const auto got = chps::nvpd_cuda(a,b).distance;
        check_cuda(cudaProfilerStop());
        if (got != expected) throw std::runtime_error("captured score differs from scalar");
        return 0;
    }
    if (std::string(argv[1]) != "model") throw std::invalid_argument("unknown study");
    check_model();
    int sms, device;
    check_cuda(cudaGetDevice(&device));
    check_cuda(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device));
    const chps::detail::MatchIndex index(a,b);
    const int stride = b.size()+1;
    int reach = 1;
    for (std::size_t row=0; row<index.entries.size(); row+=stride)
        for (int j=1; j<stride; ++j)
            if (const int k = index.entries[row+j]) reach = std::max(reach,j-k+1);
    for (bool adaptive : {false, true}) model(argv[2],a,b,index,reach,sms,expected,adaptive);
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}
