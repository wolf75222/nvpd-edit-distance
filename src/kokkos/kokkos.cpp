#include "chps/kokkos.hpp"
#include "dp_internal.hpp"

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

Result nvpd_kokkos(std::string_view a, std::string_view b, bool trace) {
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

    std::size_t current = 0;
    for (int i = 1; i <= m; ++i) {
        const std::size_t previous = current;
        current = (trace ? i : i % 2) * stride;
        const auto symbol = matches.symbol[static_cast<unsigned char>(a[i - 1])];
        const auto offset = static_cast<std::size_t>(symbol) * stride;
        // The same execution-space instance orders successive row kernels.
        Kokkos::parallel_for("NvPD row", row, KOKKOS_LAMBDA(int j) {
            if (!j) { table(current) = i; return; }
            const int k = indices(offset + j);
            const int upper = table(previous + j), diagonal = table(previous + j - 1);
            const int edit = 1 + (upper < diagonal ? upper : diagonal);
            const int match = table(previous + (k > 0 ? k - 1 : 0)) + (j - k);
            table(current + j) = edit < match ? edit : match;
        });
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

} // namespace chps
