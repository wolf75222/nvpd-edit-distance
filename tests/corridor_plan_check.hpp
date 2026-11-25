#pragma once
#include "../src/cuda/corridor_plan.hpp"
#include "chps/baselines.hpp"
#include <random>
#include <stdexcept>
#include <string>

bool check_other_corridor_geometry(std::string_view a, std::string_view b);

inline void check_corridor_plan() {
    using namespace chps::detail::corridor;
    std::mt19937 rng(14961);
    const auto dna = [&](std::size_t size) {
        std::string s(size,'A');
        for (char& c : s) c = "ACGT"[rng()%4];
        return s;
    };
    for (int size : {0,1,5,6,7,31,32,33,63,64,65,127,128,129}) {
        for (int t=0; t<48; ++t) {
            const auto a = dna(size), b = dna(size+t%9);
            const int d = chps::edlib_distance(a,b).distance;
            const int u = path_cost(a,b,static_cast<int>(a.size()+b.size()));
            const int stitched = anchored_cost(a,b,static_cast<int>(a.size()+b.size()));
            if (size<=64 && fragment_score(a,b)!=d)
                throw std::runtime_error("bit-vector fragment differs from Edlib");
            if (u<d || stitched<d || u>static_cast<int>(a.size()+b.size()))
                throw std::runtime_error("invalid constructive edit bound");
            const auto index = chunk_index(b);
            std::array<int,4> previous{};
            for (int j=0; j<=static_cast<int>(b.size()); ++j) {
                if (j) previous[dna_id(b[j-1])] = j;
                for (int symbol=0; symbol<4; ++symbol)
                    if (last(index[j/32],j,symbol)!=previous[symbol])
                        throw std::runtime_error("invalid predecessor at chunk boundary");
            }
        }
    }
    auto a = dna(33027);
    auto b = a;
    b.insert(511,"GCAA");
    b.erase(15031,3);
    b.replace(26430,8,"ACCTACCT");
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    auto* volatile build = &make_layout;
    const auto plan = build(a,b);
    if (plan.edges.size() != (b.size()/columns+1)*rows || !check_other_corridor_geometry(a,b))
        throw std::runtime_error("linked translation units mixed two CUDA tile geometries");
    if (plan.path!=Path::tiles || plan.upper<chps::edlib_distance(a,b).distance)
        throw std::runtime_error("expected a bounded NvPD corridor");
    for (int i=0; i<=m; ++i) {
        const auto [lo,hi] = output_slice(i,m,n,plan.upper);
        const int margin = (plan.upper-(n-m))/2;
        if (lo>std::max(0,i-margin) || hi<std::min(n,i+n-m+margin) || hi-lo+1>plan.capacity)
            throw std::runtime_error("rounded corridor excludes a feasible diagonal");
    }
    for (std::size_t bidx=0; bidx<plan.edges.size()/rows; ++bidx) {
        const int end = std::min(n,static_cast<int>(bidx+1)*columns-1);
        const int* edge = plan.edges.data()+bidx*rows;
        for (int t=1; t<rows; ++t) {
            for (int j=std::max(1,edge[t]); j<=end; ++j) {
                if (edge[t-1]>j-1) throw std::runtime_error("missing diagonal halo");
                for (int s=0; s<4; ++s) {
                    const int k=last(plan.index[j/32],j,s);
                    if (k && edge[t-1]>k-1) throw std::runtime_error("missing match halo");
                }
            }
        }
    }
    if (make_layout(a,a).path!=Path::length_difference ||
        make_layout(a,"T"+a).path!=Path::length_difference)
        throw std::runtime_error("failed to recognize the length lower bound");
    b[760]='N';
    if (make_layout(a,b).path!=Path::alphabet)
        throw std::runtime_error("ambiguity was encoded as DNA");
}
