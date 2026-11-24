#include "chps/fasta.hpp"
#include "local_tiles.cuh"
#include <cstdint>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 2) throw std::runtime_error("geometry requires one FASTA path");
    auto seq = chps::read_fasta(argv[1]);
    if (seq.size() != 2) throw std::runtime_error("expected two records");
    if (seq[0].bases.size() > seq[1].bases.size()) std::swap(seq[0], seq[1]);
    const int n = static_cast<int>(seq[1].bases.size());
    const std::size_t stride = static_cast<std::size_t>(n) + 1;
    chps::detail::MatchIndex index(seq[0].bases, seq[1].bases);
    int sms, limit;
    if (cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0) != cudaSuccess ||
        cudaDeviceGetAttribute(&limit, cudaDevAttrMaxSharedMemoryPerBlock, 0) != cudaSuccess)
        throw std::runtime_error("CUDA attributes unavailable");
    const int columns = stride / 1024 >= 2 * sms ? 1024 : 256;
    auto plan = chps::detail::make_local_plan(index, n, columns, limit / (2 * sizeof(int)));
    std::size_t sentinels = 0;
    int reach = 1;
    for (std::size_t row = 0; row < index.entries.size(); row += stride)
        for (int j = 0; j <= n; ++j) {
            const int last = index.entries[row + j];
            sentinels += !last || j - last >= 255;
            if (last) reach = std::max(reach, j - last + 1);
        }
    std::string path = "single-row";
    if (seq[0].bases.size() >= chps::detail::local_tile_rows) {
        if (reach <= 64) path = "uniform-tile";
        else if (!plan.left.empty())
            path = index.entries.size() == 4 * stride && 3 * sizeof(int) * plan.pitch <= limit
                ? "packed-tile" : "direct-tile";
    }
    std::cout << sms << ',' << limit << ',' << columns << ',' << plan.pitch << ',' << reach << ','
              << index.entries.size() << ',' << sentinels << ',' << path << '\n';
}
