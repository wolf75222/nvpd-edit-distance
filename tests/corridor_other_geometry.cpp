// A benchmark can override the GPU tile while its C++ callers retain defaults.
#undef CHPS_NVPD_CORRIDOR_ROWS
#undef CHPS_NVPD_CORRIDOR_COLUMNS
#define CHPS_NVPD_CORRIDOR_ROWS 32
#define CHPS_NVPD_CORRIDOR_COLUMNS 64
#include "../src/cuda/corridor_plan.hpp"

bool check_other_corridor_geometry(std::string_view a, std::string_view b) {
    using namespace chps::detail::corridor;
    // Taking the address exercises the linked symbol even in optimized builds.
    auto* volatile build = &make_layout;
    auto* volatile slice = &output_slice;
    const auto layout = build(a,b);
    const auto [lo,hi] = slice(193,static_cast<int>(a.size()),static_cast<int>(b.size()),layout.upper);
    return layout.path == Path::tiles && layout.edges.size() == (b.size()/64+1)*32
           && layout.edges[2*32-1] == 64 && lo == 128 && hi == 255;
}
