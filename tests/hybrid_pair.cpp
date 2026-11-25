#include "chps/cuda.hpp"
#include <stdexcept>
#include <string>

namespace {
using Backend = chps::HybridPairBackend;

void check(std::string_view a, std::string_view b, Backend expected) {
    if (chps::hybrid_pair_backend(a, b) != expected ||
        chps::hybrid_pair_backend(b, a) != expected)
        throw std::runtime_error("unexpected hybrid pair route");
}
}

int main(int argc, char** argv) {
    const std::string reference(100, 'A');
    const std::string six = std::string(94, 'A') + "NCCCCC";
    check(reference, six, Backend::wfa2);
    check(reference, std::string(93, 'A') + "NNCCCCC", Backend::nvpd_cuda);
    check(reference, std::string(93, 'A'), Backend::nvpd_cuda);
    check(reference, std::string(94, 'A'), Backend::wfa2);
    check("", "", Backend::wfa2);
    for (const auto length : {131072, 131073}) {
        const std::string sequence(length, 'N');
        check(sequence, sequence, length == 131072 ? Backend::wfa2 : Backend::nvpd_cuda);
    }
    Backend selected = Backend::nvpd_cuda;
    const auto score = chps::nvpd_hybrid_score(reference, six, false, &selected);
    if (selected != Backend::wfa2 || score.distance != 6 || score.device_bytes != 0)
        throw std::runtime_error("invalid CPU result");
    // Two nearby indels defeat the greedy path. The anchored fragment costs 2/40.
    const std::string_view source = "TCGTAGACAGAAATCGCGAGGAGACTAGCCAAGTTACGTA";
    const std::string_view target = "TCGTAGACAGAAATCGCGAGGGACGTAGCCAAGTTACGTA";
    check(source, target, Backend::wfa2);
    if (chps::nvpd_hybrid_score(source, target).distance != 2)
        throw std::runtime_error("invalid anchored CPU result");
    std::string query(65536, 'A');
    for (std::size_t i = 0; i < query.size(); ++i) query[i] = "ACGT"[i % 4];
    std::string changed = query;
    // Every substitution adds one C; this count difference is also a lower bound.
    int edits = 0;
    for (std::size_t i = 0; i < changed.size(); i += 10) { changed[i] = 'C'; ++edits; }
    check(query, changed, Backend::nvpd_cuda_persistent);
    check(std::string_view(query).substr(1), std::string_view(changed).substr(1), Backend::nvpd_cuda);
    std::string ambiguous = changed;
    ambiguous[1] = 'N';
    check(query, ambiguous, Backend::nvpd_cuda);
    if (argc == 2 && std::string_view(argv[1]) == "--cuda") {
        const auto gpu = chps::nvpd_hybrid_score(reference, std::string(100, 'C'), false, &selected);
        if (selected != Backend::nvpd_cuda || gpu.distance != 100)
            throw std::runtime_error("invalid CUDA result");
        const auto persistent = chps::nvpd_hybrid_score(query, changed, false, &selected);
        if (selected != Backend::nvpd_cuda_persistent || persistent.distance != edits)
            throw std::runtime_error("invalid persistent CUDA result");
    }
    try { chps::nvpd_hybrid_score(reference, six, true); }
    catch (const std::invalid_argument&) { return 0; }
    throw std::runtime_error("traceback accepted by score-only policy");
}
