#include "scheduled_matrix.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

// The final distance matrix alone cannot detect a wrong traversal or partition.
int main() {
    for (int side : {1,2,4,8,16,128}) {
        std::vector<std::pair<int,int>> order(side*side,{-1,-1});
        for (int i=0; i<side; ++i) for (int j=0; j<side; ++j) {
            auto key=chps::mpi::hilbert_index(side,i,j);
            if (key>=order.size() || order[key].first>=0) throw std::logic_error("repeated Hilbert key");
            order[key]={i,j};
        }
        for (std::size_t k=1; k<order.size(); ++k) {
            const auto [x,y]=order[k-1];
            const auto [u,v]=order[k];
            if (std::abs(x-u)+std::abs(y-v)!=1) throw std::logic_error("non-neighboring grid points");
        }
    }
    const std::vector<chps::Record> records{{"a","A"},{"b","AA"},{"c","AAA"},{"d","AAAA"}};
    std::vector<chps::mpi::Task> tasks;
    for (int i=0;i<4;++i) for (int j=i+1;j<4;++j) tasks.push_back({i,j,0,0});
    // Six pair weights [2,3,4,6,8,12]. Target prefixes for three ranks are 11 and 23.
    if (chps::mpi::weighted_cuts(tasks,records,3)!=std::vector<long long>{0,3,5,6})
        throw std::logic_error("incorrect weighted work partition");
    const std::vector<chps::Record> empty(4);
    if (chps::mpi::weighted_cuts(tasks,empty,4)!=std::vector<long long>{0,1,3,4,6})
        throw std::logic_error("incorrect empty-work partition");
    std::cout << "Hilbert order and cut positions verified\n";
}
