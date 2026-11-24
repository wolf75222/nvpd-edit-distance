#include "chps/baselines.hpp"
#include "chps/banded.hpp"
#include "chps/adaptive.hpp"
#include "chps/split.hpp"
#ifdef CHPS_CUDA
#include "chps/cuda.hpp"
#endif
#ifdef CHPS_KOKKOS
#include "chps/kokkos.hpp"
#include <optional>
#endif
#include "chps/fasta.hpp"
#include <omp.h>
#include <sys/resource.h>
#include <algorithm>
#include <charconv>
#include <chrono>
#include <climits>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using namespace chps;
using View = std::string_view;
using Kernel = Result (*)(View, View, bool);
constexpr View help =
    "editdist {pair|matrix} --algorithm {seq|hirschberg|diagonal|nvpd-omp|nvpd-adaptive|nvpd-split|nvpd-cuda|nvpd-cuda-local|nvpd-kokkos|banded|banded-trim|edlib|wfa2|biwfa}\n"
    "         --input FILE.fa [--fai FILE.fa.fai] [--threads N] [--trace]\n"
    "pair: two FASTA records, JSON output. matrix: square labeled TSV.\n"
    "--trace: pair with edlib, wfa2, biwfa, nvpd-split, hirschberg (linear space), or seq/NvPD (lengths <= 20000).\n";

struct Options {
    std::string mode, algorithm, path, fai;
    int threads = 1;
    bool trace = false;
};
Options parse(int argc, char** argv) {
    if (argc < 2) throw std::invalid_argument("missing mode");
    Options o;
    o.mode = argv[1];
    if (o.mode != "pair" && o.mode != "matrix") throw std::invalid_argument("unknown mode");
    for (int i = 2; i < argc; ++i) {
        View arg = argv[i];
        if (arg == "--trace") { o.trace = true; continue; }
        if (++i == argc) throw std::invalid_argument("missing value for " + std::string(arg));
        View value = argv[i];
        if (arg == "--algorithm") o.algorithm = value;
        else if (arg == "--input") o.path = value;
        else if (arg == "--fai") {
            if (value.empty()) throw std::invalid_argument("--fai requires an index path");
            o.fai = value;
        }
        else if (arg == "--threads") {
            const auto [end, ec] = std::from_chars(value.data(), value.data() + value.size(), o.threads);
            if (ec != std::errc{} || end != value.data() + value.size() || o.threads < 1)
                throw std::invalid_argument("--threads requires a positive integer");
        } else throw std::invalid_argument("unknown option " + std::string(arg));
    }
    if (o.path.empty() || o.algorithm.empty()) throw std::invalid_argument("input and algorithm required");
    const bool supports_trace = o.algorithm == "seq" || o.algorithm == "hirschberg" ||
                                o.algorithm == "nvpd-omp" || o.algorithm == "nvpd-cuda" || o.algorithm == "nvpd-cuda-local" ||
                                o.algorithm == "nvpd-kokkos" || o.algorithm == "nvpd-split" || o.algorithm == "edlib" ||
                                o.algorithm == "wfa2" || o.algorithm == "biwfa";
    if (o.trace && (o.mode != "pair" || !supports_trace))
        throw std::invalid_argument("unsupported traceback mode");
    return o;
}

Kernel resolve(View name) {
    if (name == "seq") return sequential;
    if (name == "hirschberg") return hirschberg;
    if (name == "nvpd-adaptive") return nvpd_adaptive;
    if (name == "nvpd-split") return nvpd_split;
    if (name == "nvpd-omp") return nvpd_omp;
    if (name == "banded" || name == "banded-trim")
        return name == "banded" ? banded : banded_trim;
    if (name == "nvpd-kokkos") {
#ifdef CHPS_KOKKOS
        return nvpd_kokkos;
#else
        throw std::invalid_argument("Kokkos backend requires -DCHPS_KOKKOS=ON");
#endif
    }
    if (name == "diagonal") return [](View a, View b, bool) { return antidiagonal(a, b); };
    if (name == "edlib") return [](View a, View b, bool trace) { return trace ? edlib_trace(a, b) : edlib_distance(a, b); };
    if (name == "wfa2") return [](View a, View b, bool trace) { return trace ? wfa_trace(a, b) : wfa_distance(a, b); };
    if (name == "biwfa") return biwfa;
    if (name == "nvpd-cuda" || name == "nvpd-cuda-local") {
#ifdef CHPS_CUDA
        return name == "nvpd-cuda-local" ? nvpd_cuda_local : nvpd_cuda;
#else
        throw std::invalid_argument("CUDA backend requires a build with -DCHPS_CUDA=ON");
#endif
    }
    throw std::invalid_argument("unknown algorithm: " + std::string(name));
}

long peak_rss_kib() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
    return usage.ru_maxrss / 1024;
#else
    return usage.ru_maxrss;
#endif
}

void run_pair(const Options& o, Kernel kernel, const std::vector<Record>& records) {
    if (records.size() != 2) throw std::invalid_argument("pair requires exactly two FASTA records");
    const auto& a = records[0].bases;
    const auto& b = records[1].bases;
    const auto start = std::chrono::steady_clock::now();
    const auto result = kernel(a, b, o.trace);
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    std::cout << std::setprecision(12) << "{\"algorithm\":\"" << o.algorithm
              << "\",\"distance\":" << result.distance << ",\"m\":" << a.size() << ",\"n\":" << b.size()
              << ",\"seconds\":" << elapsed.count() << ",\"peak_rss_kib\":" << peak_rss_kib()
              << ",\"device_bytes\":" << result.device_bytes;
#ifdef CHPS_KOKKOS
    if (o.algorithm == "nvpd-kokkos")
        std::cout << ",\"execution_space\":\"" << kokkos_execution_space() << '\"';
#endif
    if (o.trace) std::cout << ",\"cigar\":\"" << result.cigar << '"';
    std::cout << "}\n";
}

template <typename Score>
void run_matrix(const Options& o, Score&& kernel, const std::vector<Record>& records) {
    const auto n = records.size();
    std::vector<std::vector<int>> distances(n, std::vector<int>(n));
    const auto team = static_cast<std::size_t>(o.threads);
    const bool grouped = o.algorithm == "nvpd-omp" || o.algorithm == "nvpd-adaptive" || o.algorithm == "nvpd-split";
    const bool pair_parallel = grouped && team > 1 && n * (n - 1) / 2 >= team;
    const bool across = pair_parallel || (team > 1 &&
        (o.algorithm == "edlib" || o.algorithm == "wfa2" || o.algorithm == "biwfa" || o.algorithm == "seq" ||
         o.algorithm == "banded" || o.algorithm == "banded-trim" || o.algorithm == "hirschberg"));
    // Split long triangular rows into independent groups. This also exposes
    // enough work when there are fewer rows than workers but enough pairs.
    const std::size_t groups = pair_parallel ? std::min(team, n - 1) : 1;
    std::exception_ptr failure;
    auto row = [&](std::size_t i, std::size_t group) {
        for (std::size_t j = i + 1 + group; j < n; j += groups)
            distances[j][i] = distances[i][j] = kernel(records[i].bases, records[j].bases, false).distance;
    };
    if (across) {
        #pragma omp parallel
        {
            // This setting belongs to the current implicit task. Its NvPD call
            // gets one SIMD worker, including when nested teams are enabled.
            if (pair_parallel) omp_set_num_threads(1);
            #pragma omp for collapse(2) schedule(dynamic, 1)
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t group = 0; group < groups; ++group) {
                    try { row(i, group); } catch (...) {
                        #pragma omp critical
                        { if (!failure) failure = std::current_exception(); }
                    }
                }
            }
        }
        if (failure) std::rethrow_exception(failure);
    } else {
        for (std::size_t i = 0; i < n; ++i) row(i, 0);
    }
    std::cout << "id";
    for (const auto& r : records) std::cout << '\t' << r.id;
    std::cout << '\n';
    for (std::size_t i = 0; i < n; ++i) {
        std::cout << records[i].id;
        for (int distance : distances[i]) std::cout << '\t' << distance;
        std::cout << '\n';
    }
}
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && View(argv[1]) == "--help") { std::cout << help; return 0; }
        const auto o = parse(argc, argv);
        const auto kernel = resolve(o.algorithm);
        const auto records = o.fai.empty() ? chps::read_fasta(o.path) : chps::read_fasta(o.path, o.fai);
        for (const auto& r : records)
            if (r.bases.size() > (INT_MAX - 2) / 2)
                throw std::length_error("sequence exceeds the 32-bit alignment index range");
        omp_set_dynamic(0);
        omp_set_num_threads(o.threads);
#ifdef CHPS_KOKKOS
        std::optional<chps::KokkosSession> session;
        if (o.algorithm == "nvpd-kokkos") session.emplace(o.threads);
#endif
        if (o.mode == "pair") {
            run_pair(o, kernel, records);
#ifdef CHPS_CUDA
        } else if (o.algorithm == "nvpd-cuda" || o.algorithm == "nvpd-cuda-local") {
            chps::CudaWorkspace workspace(o.algorithm == "nvpd-cuda" ? chps::CudaHalo::global : chps::CudaHalo::local);
            run_matrix(o, [&](View a, View b, bool trace) { return workspace.distance(a, b, trace); }, records);
#endif
        } else {
            run_matrix(o, kernel, records);
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "editdist: " << e.what() << '\n';
        return 1;
    }
}
