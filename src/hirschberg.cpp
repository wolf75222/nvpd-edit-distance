#include "dp_internal.hpp"

#include <algorithm>
#include <numeric>
#include <utility>

namespace chps {
namespace {

using View = std::string_view;

class Path {
public:
    explicit Path(bool swapped) : insertion(swapped ? 'D' : 'I'), deletion(swapped ? 'I' : 'D') {}

    void add(char op, std::size_t length) {
        if (!length) return;
        if (op != last) flush();
        last = op;
        count += length;
        if (op != '=') cost += static_cast<int>(length);
    }

    Result finish() {
        flush();
        return {cost, std::move(cigar), 0};
    }

    const char insertion, deletion;

private:
    void flush() {
        if (!count) return;
        cigar += std::to_string(count);
        cigar += last;
        count = 0;
    }

    std::string cigar;
    char last = 0;
    std::size_t count = 0;
    int cost = 0;
};

template<class RowIterator, class ColumnIterator>
void scores(RowIterator begin, RowIterator end, ColumnIterator columns, std::size_t n, int* row) {
    std::iota(row, row + n + 1, 0);
    for (auto base = begin; base != end; ++base) {
        int diagonal = row[0];
        ++row[0];
        for (std::size_t j = 1; j <= n; ++j) {
            const int upper = row[j];
            row[j] = std::min({upper + 1, row[j-1] + 1, diagonal + (*base != columns[j-1])});
            diagonal = upper;
        }
    }
}

void one_letter(char letter, View text, char gap, Path& path) {
    const auto match = text.find(letter);
    const auto pivot = match == View::npos ? 0 : match;
    path.add(gap, pivot);
    path.add(match == View::npos ? 'X' : '=', 1);
    path.add(gap, text.size() - pivot - 1);
}

void align(View a, View b, int* forward, int* reverse, Path& path) {
    if (a.empty() || b.empty()) {
        path.add(path.insertion, b.size());
        path.add(path.deletion, a.size());
        return;
    }
    if (a.size() == 1) {
        one_letter(a.front(), b, path.insertion, path);
        return;
    }
    if (b.size() == 1) {
        one_letter(b.front(), a, path.deletion, path);
        return;
    }

    const auto middle = a.size() / 2;
    const View left = a.substr(0, middle), right = a.substr(middle);
    scores(left.begin(), left.end(), b.begin(), b.size(), forward);
    scores(right.rbegin(), right.rend(), b.rbegin(), b.size(), reverse);

    std::size_t split = 0;
    int best = forward[0] + reverse[b.size()];
    for (std::size_t j = 1; j <= b.size(); ++j) {
        const int cost = forward[j] + reverse[b.size()-j];
        if (cost < best) {
            split = j;
            best = cost;
        }
    }
    // The split is sufficient for both children; they share the two score rows.
    align(left, b.substr(0, split), forward, reverse, path);
    align(right, b.substr(split), forward, reverse, path);
}
} // namespace

Result hirschberg(View a, View b, bool trace) {
    if (!trace) return sequential(a, b, false);
    detail::validate(a, b, false);
    if (a.empty() || b.empty()) return detail::empty_result(a, b, true);
    const bool swapped = a.size() < b.size();
    if (swapped) std::swap(a, b);
    Path path(swapped);
    const auto stride = b.size() + 1;
    std::vector<int> rows(2 * stride);
    align(a, b, rows.data(), rows.data() + stride, path);
    return path.finish();
}
} // namespace chps
