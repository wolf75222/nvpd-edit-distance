#include "chps/baselines.hpp"
#include "chps/split.hpp"
#include "cigar_replay.hpp"
#include <omp.h>
#include <random>

chps::SplitStats check(std::string_view a, std::string_view b) {
    chps::SplitStats stats;
    auto result = chps::nvpd_split_profile(a, b, true, stats);
    replay(a, b, result);
    if (result.distance != chps::edlib_distance(a, b).distance ||
        result.distance != chps::nvpd_split(a, b, false).distance)
        throw std::runtime_error("split distance differs from reference");
    return stats;
}

int main() {
    omp_set_num_threads(2);
    check("", ""); check("A", ""); check("NRY", "YN"); check("CAT", "ATC");
    std::mt19937 rng(42);
    std::string a(131072, 'A');
    for (auto& c : a) c = "ACGT"[rng() % 4];
    std::string b = a;
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
