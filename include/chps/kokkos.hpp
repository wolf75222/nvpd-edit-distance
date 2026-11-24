#pragma once

#include "chps/dp.hpp"

namespace chps {

// Create before calling the backend. An embedding application may own Kokkos instead.
class KokkosSession {
public:
    explicit KokkosSession(int threads);
    ~KokkosSession();
    KokkosSession(const KokkosSession&) = delete;
    KokkosSession& operator=(const KokkosSession&) = delete;

private:
    bool owner_;
};

const char* kokkos_execution_space();
Result nvpd_kokkos(std::string_view a, std::string_view b, bool trace = false);

} // namespace chps
