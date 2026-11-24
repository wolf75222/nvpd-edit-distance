#include "chps/baselines.hpp"
#include "../src/cpu/sparse_frontier.hpp"
#include <algorithm>
#include <vector>
#include "chps/split.hpp"
#include "cigar_replay.hpp"
#include <omp.h>
#include <random>

chps::SplitStats check(std::string_view a, std::string_view b) {
    chps::SplitStats stats;
    auto result = chps::nvpd_split_profile(a, b, true, stats);
    replay(a, b, result);
    if (result.distance != chps::edlib_distance(a, b).distance ||
        result.distance != chps::nvpd_split(a, b, false).distance ||
        result.distance != chps::nvpd_sparse(a, b, false).distance)
        throw std::runtime_error("split distance differs from reference");
    return stats;
}

chps::SplitStats check_filter(std::string_view a, std::string_view b) {
    const chps::SplitComposition first(a), second(b);
    const auto lower = first.lower_bound(second);
    const auto oracle = chps::edlib_distance(a,b).distance;
    if (lower != second.lower_bound(first) || lower > static_cast<std::size_t>(oracle))
        throw std::runtime_error("invalid composition lower bound");
    chps::SplitStats stats;
    for (bool trace : {false,true}) {
        auto path = chps::nvpd_split_profile(a,b,trace,stats,chps::SplitPolicy::retain_verified,lower);
        if (path.distance != oracle) throw std::runtime_error("filtered split distance differs from edlib");
        if (trace) replay(a,b,path);
    }
    return stats;
}

int main() {
    std::vector<std::string> words{""};
    for (std::size_t k = 0; k < words.size(); ++k) {
        if (words[k].size() == 5) break;
        for (char x : std::string_view("AN")) words.push_back(words[k] + x);
    }
    for (const auto& a : words) for (const auto& b : words) {
        const int d = chps::edlib_distance(a,b).distance;
        for (int limit : {0, std::max(0,d-1), d, d+1})
            if (chps::detail::sparse_frontier(a,b,limit) != (d <= limit ? d : -1))
                throw std::runtime_error("bounded frontier differs from the global oracle");
    }
    for (int edits : {127,128,129}) {
        const std::string a(1000,'A');
        auto b = a;
        b.replace(100, edits, std::string(edits,'N'));
        if (chps::detail::sparse_frontier(a,b,128) != (edits <= 128 ? edits : -1))
            throw std::runtime_error("sparse budget boundary differs from expected result");
    }
    omp_set_num_threads(2);
    check("", ""); check("A", ""); check("NRY", "YN"); check("CAT", "ATC");
    std::mt19937 rng(42);
    std::string a(131072, 'A');
    for (auto& c : a) c = "ACGT"[rng() % 4];
    std::string b = a;
    const auto small = a.substr(0,16384);
    for (int edits : {127,128,1025}) {
        auto mutant = small;
        mutant.replace(0,edits,std::string(edits,'N'));
        const auto stats = check_filter(small,mutant);
        if (stats.composition_pruned != (edits >= 128))
            throw std::runtime_error("incorrect composition filter at certificate boundary");
    }
    const auto longer = std::string(129,'N') + small;
    for (bool reverse : {false,true}) {
        const auto stats = reverse ? check_filter(longer,small) : check_filter(small,longer);
        if (stats.composition_pruned || stats.accepted_cuts != 1)
            throw std::runtime_error("pruned a cut certified by the length difference");
    }
    b.replace(31007, 1, "NN"); b.erase(71907, 3); b.insert(98007, "Y");
    if (check(a, b).accepted_cuts != 7 || check(b, a).accepted_cuts != 7)
        throw std::runtime_error("missing exact split in unique sequence");
    b = a; b.replace(4000, 1300, std::string(1300, 'N'));
    if (!check(a, b).global_retry) throw std::runtime_error("missing global retry");
    // Six cuts remain usable when a local difference interrupts one witness interval.
    b = a;
    b[16384-80] = 'N';
    for (int pos = 300; pos < 390; pos += 10) b[pos] = 'N';
    auto partial = check(a, b);
    if (partial.accepted_cuts != 6 || partial.proposed_cuts != 7 || partial.global_retry)
        throw std::runtime_error("independent certification lost usable cuts");
    // The only matching seed in B is a displaced duplicate; forcing it costs at least 512.
    auto duplicated = a.substr(0, 16384);
    duplicated.replace(8192+256-64, 128, duplicated.substr(8192-64, 128));
    auto changed = duplicated;
    changed[8192] = 'N';
    auto false_match = check(duplicated, changed);
    if (false_match.proposed_cuts != 1 || false_match.accepted_cuts != 0 || !false_match.global_retry ||
        chps::nvpd_split(duplicated, changed).distance != 1)
        throw std::runtime_error("displaced duplicate forced a nonoptimal path");
    std::string scarce(16384, 'A');
    scarce.replace(8192-64, 128, a.substr(8192-64,128));
    auto eight_edits = scarce;
    for (int i=0; i<8; ++i) eight_edits[i*100] = 'N';
    auto repeated_words = check(scarce, eight_edits);
    if (repeated_words.proposed_cuts != 1 || !repeated_words.global_retry)
        throw std::runtime_error("failed to reject nonunique witness words");
    const std::string repeat(66000, 'A');
    if (check(repeat, repeat+"T").proposed_cuts) throw std::runtime_error("nonunique seed");
}
