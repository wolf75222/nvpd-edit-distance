// Export NvPD's actual corridor plan without launching or timing a solver.
#include "../../../src/cuda/corridor_plan.hpp"
#include "chps/fasta.hpp"
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace chps::detail::corridor;

namespace {
const char* name(Path path) {
    switch (path) {
        case Path::size: return "dense-size";
        case Path::upper_bound: return "dense-cost-limit";
        case Path::length_difference: return "length-difference";
        case Path::alphabet: return "dense-alphabet";
        case Path::wide_halo: return "dense-halo";
        case Path::tiles: return "corridor";
    }
    throw std::logic_error("unknown corridor path");
}

std::uint64_t cells(const Layout& layout, int m, int n) {
    const int blocks = n/columns+1;
    std::vector<std::uint64_t> prefix(static_cast<std::size_t>(blocks)+1);
    for (int block=0; block<blocks; ++block) {
        const int end = std::min(n,(block+1)*columns-1);
        prefix[block+1] = prefix[block];
        for (int t=0; t<rows; ++t)
            prefix[block+1] += end-layout.edges[static_cast<std::size_t>(block)*rows+t]+1;
    }
    std::uint64_t count = 0;
    for (int i=1; i<=m; i+=rows) {
        const int depth = std::min(rows,m-i+1);
        const auto slice = output_slice(i+depth-1,m,n,layout.upper);
        const int first = slice.first/columns, last_block = slice.last/columns;
        if (depth==rows) {
            count += prefix[last_block+1]-prefix[first];
        } else {
            for (int block=first; block<=last_block; ++block)
                for (int t=rows-depth; t<rows; ++t)
                    count += std::min(n,(block+1)*columns-1)
                        -layout.edges[static_cast<std::size_t>(block)*rows+t]+1;
        }
    }
    return count;
}
}

int main(int argc, char** argv) {
    if (argc<2) {std::cerr<<"usage: corridor_plan pair.fa [...]\n";return 2;}
    std::cout<<"path,rows,columns,policy,upper,frontier_capacity,shared_pitch,launches,tile_cells,device_bytes\n";
    for (int i=1; i<argc; ++i) {
        const auto records = chps::read_fasta(argv[i]);
        if (records.size()!=2) throw std::invalid_argument("expected two FASTA records");
        std::string_view a=records[0].bases,b=records[1].bases;
        if (a.size()>b.size()) std::swap(a,b);
        if (a.size()+b.size()>=static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("input exceeds index range");
        const int m=static_cast<int>(a.size()), n=static_cast<int>(b.size());
        const auto layout=make_layout(a,b);
        std::cout<<std::quoted(std::filesystem::path(argv[i]).filename().string(),'"','"')
                 <<','<<m<<','<<n<<','<<name(layout.path)<<',';
        if (layout.upper>=0 && layout.path!=Path::upper_bound) std::cout<<layout.upper;
        if (layout.path==Path::tiles) {
            const auto bytes=layout.index.size()*sizeof(Chunk)
                +(layout.edges.size()+2u*layout.capacity)*sizeof(int);
            std::cout<<','<<layout.capacity<<','<<layout.pitch<<','<<(m+rows-1)/rows
                     <<','<<cells(layout,m,n)<<','<<bytes;
        } else std::cout<<",,,,,";
        std::cout<<'\n';
    }
}
