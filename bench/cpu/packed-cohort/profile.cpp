#include "chps/baselines.hpp"
#include "chps/fasta.hpp"
#include "chps/packed.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 5) throw std::invalid_argument("profile FASTA METHOD PERF_CONTROL_FIFO PERF_ACK_FIFO");
    auto rows = chps::read_fasta(argv[1]);
    if (rows.size()<128) throw std::invalid_argument("expected at least 128 records");
    rows.resize(128);
    std::vector<chps::PackedSequence> compact;
    for (const auto& row : rows) compact.emplace_back(row.bases);
    const std::string method = argv[2];
    auto kernel = [&](int i,int j) {
        if (method=="byte") return chps::byte_wavefront_score(rows[i].bases, rows[j].bases);
        if (method=="packed") return chps::packed_score(compact[i],compact[j]);
        if (method=="one-sided") return chps::packed_score(compact[i],compact[j],chps::PackedFront::one_sided);
        if (method=="wfa2") return chps::wfa_distance(rows[i].bases,rows[j].bases).distance;
        if (method=="edlib") return chps::edlib_distance(rows[i].bases,rows[j].bases).distance;
        throw std::invalid_argument("unknown solver " + method);
    };
    auto matrix = [&] {
        std::vector<int> scores;
        scores.reserve(rows.size()*(rows.size()-1)/2);
        for (int i=0; i<128; ++i)
            for (int j=i+1; j<128; ++j) scores.push_back(kernel(i,j));
        return scores;
    };
    auto warm = matrix();
    std::ofstream ctl(argv[3]); std::ifstream ack(argv[4]);
    auto signal = [&](const char* command) {
        ctl << command << std::endl;
        std::string answer;
        if (!std::getline(ack,answer)) throw std::runtime_error("missing perf acknowledgement");
        answer.erase(std::remove(answer.begin(),answer.end(),'\0'),answer.end());
        if (answer != "ack") throw std::runtime_error("unexpected perf acknowledgement: " + answer);
    };
    signal("enable");
    const auto before=std::chrono::steady_clock::now();
    auto measured=matrix();
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
    signal("disable");
    if (warm!=measured) throw std::runtime_error("warm and timed matrices disagree");
    std::uint64_t hash=UINT64_C(14695981039346656037);
    for (const auto value:measured) hash=(hash ^ static_cast<std::uint64_t>(value))*UINT64_C(1099511628211);
    std::cout << hash << ',' << elapsed << '\n';
}
