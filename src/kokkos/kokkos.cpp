#include "chps/kokkos.hpp"
#include "dp_internal.hpp"
#include "row_policy.hpp"

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace chps {
namespace {
using Exec = Kokkos::DefaultExecutionSpace;
using Memory = Exec::memory_space;
using Index = Kokkos::RangePolicy<Exec, Kokkos::IndexType<int>>;
using Buffer = Kokkos::View<int*, Memory>;

// One host team stays in the same parallel region for all rows. Each worker
// owns contiguous columns, and the team barrier makes the next row readable.
template<class HostExec>
bool host_team_rows(const HostExec& exec, const char* sequence, const int* symbols,
                    const int* indices, int* table, int m, int n, bool trace) {
    using Policy = Kokkos::TeamPolicy<HostExec>;
    using Member = typename Policy::member_type;
    const std::size_t stride = static_cast<std::size_t>(n) + 1;
    const auto compute = [=](const Member& team) {
        const int first = static_cast<int>(static_cast<std::size_t>(n) * team.team_rank() / team.team_size()) + 1;
        const int end = static_cast<int>(static_cast<std::size_t>(n) * (team.team_rank() + 1) / team.team_size()) + 1;
        for (int i = 1; i <= m; ++i) {
            const int* previous = table + (trace ? i - 1 : (i - 1) % 2) * stride;
            int* current = table + (trace ? i : i % 2) * stride;
            const auto offset = symbols[static_cast<unsigned char>(sequence[i - 1])] * stride;
            if (team.team_rank() == 0) current[0] = i;
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, first, end), [&](int j) {
                const int k = indices[offset + j];
                const int upper = previous[j], diagonal = previous[j - 1];
                const int edit = 1 + (upper < diagonal ? upper : diagonal);
                const int match = previous[k > 0 ? k - 1 : 0] + (j - k);
                current[j] = edit < match ? edit : match;
            });
            team.team_barrier();
        }
    };
    const Policy automatic(exec, 1, Kokkos::AUTO);
    const int workers = exec.concurrency();
    if (workers > automatic.team_size_max(compute, Kokkos::ParallelForTag{})) return false;
    Kokkos::parallel_for("NvPD persistent host team", Policy(exec, 1, workers), compute);
    return true;
}
} // namespace

KokkosSession::KokkosSession(int threads) : owner_(!Kokkos::is_initialized()) {
    if (threads < 1) throw std::invalid_argument("Kokkos needs at least one host thread");
    if (owner_) {
        if (Kokkos::is_finalized()) throw std::logic_error("Kokkos already finalized");
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(threads);
        Kokkos::initialize(settings);
    }
}

KokkosSession::~KokkosSession() {
    if (owner_) Kokkos::finalize();
}

const char* kokkos_execution_space() {
    return Exec::name();
}

namespace {
// Evaluate an intermediate row directly from the previous completed row. A
// second-row cell needs at most three such values, regardless of match reach.
KOKKOS_INLINE_FUNCTION
int intermediate(const int* prev, const int* last, int i, int j) {
    if (j == 0) return i;
    const int k = last[j];
    const int diagonal = prev[j - 1];
    if (k == j) return diagonal;
    const int upper = prev[j];
    const int edit = 1 + (upper < diagonal ? upper : diagonal);
    if (k == 0) return edit;
    const int match = prev[k - 1] + (j - k);
    return edit < match ? edit : match;
}

std::size_t fused_rows(const Index& row, std::string_view a, const detail::MatchIndex& matches,
                       const Buffer& indices, const Buffer& table, std::size_t stride) {
    std::size_t current = 0;
    for (int i = 1; i <= static_cast<int>(a.size()); i += 2) {
        const auto* prev = table.data() + current;
        current = stride - current;
        auto* out = table.data() + current;
        const int s1 = matches.symbol[static_cast<unsigned char>(a[i - 1])];
        const int* first = indices.data() + s1 * stride;
        const bool two = i < static_cast<int>(a.size());
        const int s2 = two ? matches.symbol[static_cast<unsigned char>(a[i])] : s1;
        const int* second = indices.data() + s2 * stride;
        Kokkos::parallel_for("NvPD two rows", row, KOKKOS_LAMBDA(int j) {
            if (!two) { out[j] = intermediate(prev, first, i, j); return; }
            if (!j) { out[0] = i + 1; return; }
            const int k = second[j];
            const int diagonal = intermediate(prev, first, i, j - 1);
            if (k == j) { out[j] = diagonal; return; }
            const int upper = intermediate(prev, first, i, j);
            const int edit = 1 + (upper < diagonal ? upper : diagonal);
            if (!k) { out[j] = edit; return; }
            const int match = intermediate(prev, first, i, k - 1) + (j - k);
            out[j] = edit < match ? edit : match;
        });
    }
    return current;
}

Result compute_distance(std::string_view a, std::string_view b, bool trace, bool persistent) {
    detail::validate(a, b, trace);
    if (!Kokkos::is_initialized()) throw std::logic_error("initialize Kokkos before computing distances");
    if (a.empty() || b.empty()) return detail::empty_result(a, b, trace);
    if (!trace && a.size() > b.size()) std::swap(a, b);
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    const std::size_t stride = b.size() + 1;
    const std::size_t rows = trace ? a.size() + 1 : 2;
    const detail::MatchIndex matches(a, b);
    Exec exec;
    Buffer indices(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "match indices"),
                   matches.entries.size());
    Buffer table(Kokkos::view_alloc(exec, Kokkos::WithoutInitializing, "score rows"), rows * stride);
    using HostIndices = Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
    Kokkos::deep_copy(exec, indices, HostIndices(matches.entries.data(), matches.entries.size()));
    const Index row(exec, 0, n + 1);
    Kokkos::parallel_for("NvPD first row", row, KOKKOS_LAMBDA(int j) { table(j) = j; });

    bool complete = false;
    std::size_t current = (trace ? m : m % 2) * stride;
    if (persistent) {
        if constexpr (std::is_same_v<Memory, Kokkos::HostSpace>) {
            complete = host_team_rows(Kokkos::DefaultHostExecutionSpace{}, a.data(), matches.symbol.data(),
                                      indices.data(), table.data(), m, n, trace);
        } else if (!trace) {
            current = fused_rows(row, a, matches, indices, table, stride);
            complete = true;
        }
    }
    if (!complete) {
        for (int i = 1; i <= m; ++i) {
            const auto previous = (trace ? i - 1 : (i - 1) % 2) * stride;
            current = (trace ? i : i % 2) * stride;
            const int symbol = matches.symbol[static_cast<unsigned char>(a[i - 1])];
            const auto offset = symbol * stride;
            Kokkos::parallel_for("NvPD row", row, KOKKOS_LAMBDA(int j) {
                if (!j) { table(current) = i; return; }
                const int k = indices(offset + j);
                const int upper = table(previous + j), diagonal = table(previous + j - 1);
                const int edit = 1 + (upper < diagonal ? upper : diagonal);
                const int match = table(previous + (k > 0 ? k - 1 : 0)) + (j - k);
                table(current + j) = edit < match ? edit : match;
            });
        }
    }
    exec.fence("NvPD last row ready");
    const auto bytes = std::is_same_v<Memory, Kokkos::HostSpace> ? std::size_t{0}
                       : (indices.extent(0) + table.extent(0)) * sizeof(int);
    if (trace) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, table);
        return {host(current + n), detail::traceback(host.data(), a, b), bytes};
    }
    int distance;
    Kokkos::deep_copy(distance, Kokkos::subview(table, current + n));
    return {distance, "", bytes};
}

} // namespace

Result detail::nvpd_kokkos_rows(std::string_view a, std::string_view b, bool trace) {
    return compute_distance(a, b, trace, false);
}

Result detail::nvpd_kokkos_opt(std::string_view a, std::string_view b, bool trace) {
    return compute_distance(a, b, trace, true);
}

Result nvpd_kokkos(std::string_view a, std::string_view b, bool trace) {
    return compute_distance(a, b, trace, true);
}

} // namespace chps
