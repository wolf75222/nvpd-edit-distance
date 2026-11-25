#include "chps/split.hpp"
#include "chps/baselines.hpp"
#include "dp_internal.hpp"
#include "sparse_frontier.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <exception>
#include <omp.h>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <vector>

namespace chps {
namespace {
using View = std::string_view;
using Size = std::size_t;
constexpr Size word = 16, radius = 2048, max_bound = 1024;
struct Cut { Size i, j; };
struct Piece {
    Result result{};
    Size prefix = 0, suffix = 0;
    bool dense = false, wavefront = false;
};

// CIGAR operations at either side of a cut can belong to the same run.
struct Runs {
    std::string data;
    char last = 0;
    Size count = 0;
    void flush() {
        if (count) { data += std::to_string(count); data += last; count = 0; }
    }
    void add(char op, Size length) {
        if (!length) return;
        if (op != last) { flush(); last = op; }
        count += length;
    }
    void append(View cigar) {
        while (!cigar.empty()) {
            Size n = 0;
            auto [end, error] = std::from_chars(cigar.data(), cigar.data() + cigar.size(), n);
            if (error != std::errc{} || end == cigar.data() + cigar.size() || !n)
                throw std::logic_error("invalid component CIGAR");
            add(*end, n);
            cigar.remove_prefix(end - cigar.data() + 1);
        }
    }
};

Size common_end(View a, View b, bool reverse) {
    const Size n = std::min(a.size(), b.size());
    Size k = 0;
    constexpr Size block = 512;
    while (n - k >= block) {
        const Size x = reverse ? a.size() - k - block : k;
        const Size y = reverse ? b.size() - k - block : k;
        if (std::memcmp(a.data() + x, b.data() + y, block)) break;
        k += block;
    }
    while (k < n && a[reverse ? a.size()-1-k : k] == b[reverse ? b.size()-1-k : k]) ++k;
    return k;
}

Piece align_piece(View a, View b, bool trace, bool sparse) {
    const Size prefix = common_end(a, b, false);
    a.remove_prefix(prefix); b.remove_prefix(prefix);
    const Size suffix = common_end(a, b, true);
    a.remove_suffix(suffix); b.remove_suffix(suffix);
    Piece p;
    p.prefix = prefix;
    p.suffix = suffix;
    if (a.empty() || b.empty()) {
        p.result = detail::empty_result(a, b, trace);
    } else if (std::max(a.size(), b.size()) <= 20000 && a.size() <= 262144 / b.size()) {
        p.dense = true;
        p.result = nvpd_omp(a, b, trace);
    } else {
        const int frontier = sparse && !trace ? detail::sparse_frontier(a, b, 128) : -1;
        if (frontier >= 0) {
            p.result = {frontier, {}, 0};
        } else {
            p.wavefront = true;
            p.result = biwfa(a, b, trace);
        }
    }
    if (trace) {
        Runs runs;
        runs.add('=', prefix); runs.append(p.result.cigar); runs.add('=', suffix); runs.flush();
        p.result.cigar = std::move(runs.data);
    }
    return p;
}

std::vector<Cut> candidate_cuts(View a, View b) {
    std::vector<Cut> cuts{{0,0}};
    const Size parts = std::min<Size>(8, std::min(a.size(), b.size()) / 8192);
    for (Size k = 1; k < parts; ++k) {
        const Size i = a.size() * k / parts, predicted = b.size() * k / parts;
        const Size low = predicted > radius ? predicted - radius : 0;
        const Size high = std::min(b.size(), predicted + radius + word);
        const View region = b.substr(low, high-low), key = a.substr(i, word);
        const Size hit = region.find(key);
        if (hit == View::npos || region.find(key, hit+1) != View::npos) continue;
        const Size j = low + hit;
        if (j <= cuts.back().j || j < 64 || j + 64 > b.size()) continue;
        if (a.substr(i-64,128) == b.substr(j-64,128)) cuts.push_back({i,j});
    }
    cuts.push_back({a.size(), b.size()});
    return cuts;
}

int letter(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}

// Count exact 16-mers over the union of every target position reachable at cost U.
// Ambiguous symbols break a word; their equality is still literal in the alignment.
bool witnesses(View a, View b, Cut c, Size bound, Size span) {
    const Size start = c.i - span;
    if (a.substr(start,2*span) != b.substr(c.j-span,2*span)) return false;
    const Size low = start > bound ? start-bound : 0;
    const Size high = std::min(b.size(), c.i+span+bound);
    std::unordered_map<std::uint32_t, unsigned char> counts;
    counts.reserve(2*span/word);
    for (Size p = start; p < start + 2*span; p += word) {
        std::uint32_t key = 0;
        bool clean = true;
        for (char ch : a.substr(p,word)) {
            const int x = letter(ch);
            if (x < 0) { clean = false; break; }
            key = (key << 2) | static_cast<unsigned>(x);
        }
        if (clean) counts.emplace(key, 0);
    }
    std::uint32_t value = 0;
    Size valid = 0;
    for (Size pos = low; pos < high; ++pos) {
        const int x = letter(b[pos]);
        if (x < 0) { valid = 0; value = 0; continue; }
        value = (value << 2) | static_cast<unsigned>(x);
        if (++valid >= word) {
            auto it = counts.find(value);
            if (it != counts.end() && it->second < 2) ++it->second;
        }
    }
    for (Size side : {start, c.i}) {
        Size unique = 0;
        for (Size p = side; p < side+span && unique <= bound; p += word) {
            value = 0;
            bool clean = true;
            for (char ch : a.substr(p,word)) {
                const int x = letter(ch);
                if (x < 0) { clean = false; break; }
                value = (value << 2) | static_cast<unsigned>(x);
            }
            if (clean && counts.at(value) == 1) ++unique;
        }
        if (unique <= bound) return false;
    }
    return true;
}

Size clearance(const std::vector<Cut>& cuts, Size k) {
    const Cut prev = cuts[k-1], c = cuts[k], next = cuts[k+1];
    return std::min({c.i-prev.i, c.j-prev.j, next.i-c.i, next.j-c.j});
}

bool certified_cut(View a, View b, const std::vector<Cut>& cuts, Size k, Size bound) {
    if (bound == std::max(a.size(), b.size()) - std::min(a.size(), b.size())) return true;
    if (bound > max_bound) return false;
    const Size span = 2 * (bound + 1) * word;
    // Disjoint neighborhoods let all verified optimal subpaths be spliced simultaneously.
    return clearance(cuts, k) >= 2*span && witnesses(a, b, cuts[k], bound, span);
}

Piece align_range(View a, View b, Cut first, Cut last, bool trace, bool sparse) {
    return align_piece(a.substr(first.i, last.i-first.i), b.substr(first.j, last.j-first.j), trace, sparse);
}

// Equal ends can be removed even when the proposed cut through the middle fails.
Piece merged_score(View a, View b, const std::vector<Cut>& cuts, const std::vector<Piece>& pieces,
                   Size first, Size last, bool sparse) {
    while (first < last && pieces[first].result.distance == 0) ++first;
    while (first < last && pieces[last-1].result.distance == 0) --last;
    if (first == last) return {};
    if (first + 1 == last) {
        Piece p;
        p.result.distance = pieces[first].result.distance;
        return p;
    }
    Cut from = cuts[first], to = cuts[last];
    from.i += pieces[first].prefix;
    from.j += pieces[first].prefix;
    to.i -= pieces[last-1].suffix;
    to.j -= pieces[last-1].suffix;
    return align_range(a, b, from, to, false, sparse);
}
} // namespace

SplitComposition::SplitComposition(std::string_view sequence) {
    for (unsigned char symbol : sequence) ++counts_[symbol];
}

std::size_t SplitComposition::lower_bound(const SplitComposition& other) const {
    std::size_t surplus = 0, deficit = 0;
    for (std::size_t c = 0; c < counts_.size(); ++c) {
        if (counts_[c] > other.counts_[c]) surplus += counts_[c] - other.counts_[c];
        else deficit += other.counts_[c] - counts_[c];
    }
    return std::max(surplus, deficit);
}

static Result split_impl(View a, View b, bool trace, SplitStats& stats,
                         SplitPolicy policy, Size lower_bound, bool sparse) {
    detail::validate(a, b, false);
    stats = {};
    const Size imbalance = std::max(a.size(), b.size()) - std::min(a.size(), b.size());
    auto single = [&] {
        auto p = align_piece(a, b, trace, sparse);
        stats.dense_pieces = p.dense;
        stats.wavefront_pieces = p.wavefront;
        return std::move(p.result);
    };
    if (policy == SplitPolicy::single_piece || (policy == SplitPolicy::retain_verified && imbalance > max_bound))
        return single();
    // The feasible split score U is at least the composition bound. If L exceeds
    // the length gap, U cannot use the length-equality shortcut.
    const bool filter = policy == SplitPolicy::retain_verified && lower_bound > imbalance;
    if (filter && lower_bound > max_bound) {
        stats.composition_pruned = true;
        return single();
    }
    const auto cuts = candidate_cuts(a, b);
    stats.proposed_cuts = static_cast<int>(cuts.size()) - 2;
    if (filter && stats.proposed_cuts > 0) {
        const Size min_gap = 4 * (lower_bound + 1) * word;
        Size max_gap = 0;
        for (Size k = 1; k + 1 < cuts.size(); ++k)
            max_gap = std::max(max_gap, clearance(cuts, k));
        if (max_gap < min_gap) {
            stats.composition_pruned = true;
            return single();
        }
    }
    std::vector<Piece> pieces(cuts.size() - 1);
    std::exception_ptr failure;
    #pragma omp parallel for if(pieces.size()>1 && !omp_in_parallel()) schedule(static)
    for (int k=0; k<static_cast<int>(pieces.size()); ++k) {
        try {
            pieces[k] = align_range(a, b, cuts[k], cuts[k+1], trace, sparse);
        } catch (...) {
            #pragma omp critical(chps_split_error)
            { if (!failure) failure = std::current_exception(); }
        }
    }
    if (failure) std::rethrow_exception(failure);
    int distance = 0;
    for (const Piece& p : pieces) {
        distance += p.result.distance;
        stats.dense_pieces += p.dense;
        stats.wavefront_pieces += p.wavefront;
    }
    if (pieces.size()==1) return std::move(pieces.front().result);
    const int upper = distance;
    distance = 0;
    Size first = 0;
    Runs path;
    for (Size last = 1; last < cuts.size(); ++last) {
        if (last + 1 < cuts.size()) {
            if (!certified_cut(a, b, cuts, last, upper)) {
                if (policy == SplitPolicy::require_all) {
                    stats.global_retry = true;
                    stats.accepted_cuts = 0;
                    return biwfa(a, b, trace);
                }
                continue;
            }
            ++stats.accepted_cuts;
        }
        if (last != first + 1) {
            pieces[first] = trace ? align_range(a, b, cuts[first], cuts[last], true, sparse)
                                  : merged_score(a, b, cuts, pieces, first, last, sparse);
            stats.dense_pieces += pieces[first].dense;
            stats.wavefront_pieces += pieces[first].wavefront;
        }
        distance += pieces[first].result.distance;
        if (trace) path.append(pieces[first].result.cigar);
        first = last;
    }
    stats.global_retry = stats.accepted_cuts == 0;
    path.flush();
    return {distance,std::move(path.data),0};
}

Result nvpd_split_profile(View a, View b, bool trace, SplitStats& stats,
                          SplitPolicy policy, Size lower_bound) {
    return split_impl(a, b, trace, stats, policy, lower_bound, false);
}

Result nvpd_sparse(View a, View b, bool trace) {
    SplitStats stats;
    return split_impl(a, b, trace, stats, SplitPolicy::retain_verified, 0, true);
}

Result nvpd_split(View a, View b, bool trace) {
    SplitStats stats;
    return nvpd_split_profile(a,b,trace,stats);
}
} // namespace chps
