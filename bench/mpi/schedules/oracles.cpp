#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

// Compile with the same pinned libraries used by the command-line baselines.
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    auto data = chps::read_fasta(argv[1]);
    std::sort(data.begin(),data.end(),[](const auto& l, const auto& r){return l.id<r.id;});
    std::vector<std::pair<int,int>> jobs;
    for (int i=0; i<static_cast<int>(data.size()); ++i)
        for (int j=i+1; j<static_cast<int>(data.size()); ++j) jobs.emplace_back(i,j);
    std::vector<int> ed(jobs.size()), wf(jobs.size());
    #pragma omp parallel for schedule(dynamic,1)
    for (int k=0; k<static_cast<int>(jobs.size()); ++k) {
        const auto [i,j]=jobs[k];
        ed[k]=chps::edlib_distance(data[i].bases,data[j].bases).distance;
        wf[k]=chps::wfa_distance(data[i].bases,data[j].bases).distance;
    }
    std::cout << "id_a\tid_b\tedlib\twfa2\n";
    for (std::size_t k=0; k<jobs.size(); ++k) {
        const auto [i,j]=jobs[k];
        if (ed[k]!=wf[k]) throw std::runtime_error("independent scores disagree");
        std::cout << data[i].id << '\t' << data[j].id << '\t' << ed[k] << '\t' << wf[k] << '\n';
    }
}
