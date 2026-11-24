#include "chps/baselines.hpp"
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
    "editdist {pair|matrix} --algorithm {seq|diagonal|nvpd-omp|nvpd-cuda|edlib|wfa2}\n"
    "         --input FILE.fa [--threads N] [--trace]\n"
    "pair: two FASTA records, JSON output. matrix: square labeled TSV.\n"
    "--trace: pair mode with seq, nvpd-omp or nvpd-cuda, lengths <= 20000.\n";

struct Options {
    std::string mode, algorithm, path;
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
        else if (arg == "--threads") {
            const auto [end, ec] = std::from_chars(value.data(), value.data() + value.size(), o.threads);
            if (ec != std::errc{} || end != value.data() + value.size() || o.threads < 1)
                throw std::invalid_argument("--threads requires a positive integer");
        } else throw std::invalid_argument("unknown option " + std::string(arg));
    }
    if (o.path.empty() || o.algorithm.empty()) throw std::invalid_argument("input and algorithm required");
    if (o.trace && (o.mode != "pair" || (o.algorithm != "seq" && o.algorithm != "nvpd-omp"
                                            && o.algorithm != "nvpd-cuda")))
        throw std::invalid_argument("unsupported traceback mode");
    return o;
}

Kernel resolve(View name) {
    if (name == "seq") return sequential;
    if (name == "nvpd-omp") return nvpd_omp;
    if (name == "diagonal") return [](View a, View b, bool) { return antidiagonal(a, b); };
    if (name == "edlib") return [](View a, View b, bool) { return edlib_distance(a, b); };
    if (name == "wfa2") return [](View a, View b, bool) { return wfa_distance(a, b); };
    if (name == "nvpd-cuda") {
#ifdef CHPS_CUDA
        return nvpd_cuda;
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
    if (o.trace) std::cout << ",\"cigar\":\"" << result.cigar << '"';
    std::cout << "}\n";
}

void run_matrix(const Options& o, Kernel kernel, const std::vector<Record>& records) {
    const auto n = records.size();
    std::vector<std::vector<int>> distances(n, std::vector<int>(n));
    // Serial kernels parallelize over independent rows. OpenMP kernels own their team.
    const bool across = o.threads > 1 && (o.algorithm == "edlib" || o.algorithm == "wfa2" || o.algorithm == "seq");
    std::exception_ptr failure;
    auto row = [&](std::size_t i) {
        for (std::size_t j = i + 1; j < n; ++j)
            distances[j][i] = distances[i][j] = kernel(records[i].bases, records[j].bases, false).distance;
    };
    if (across) {
#pragma omp parallel for schedule(dynamic, 1)
        for (std::size_t i = 0; i < n; ++i) {
            try { row(i); } catch (...) {
#pragma omp critical
                { if (!failure) failure = std::current_exception(); }
            }
        }
        if (failure) std::rethrow_exception(failure);
    } else {
        for (std::size_t i = 0; i < n; ++i) row(i);
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
        const auto records = chps::read_fasta(o.path);
        for (const auto& r : records)
            if (r.bases.size() > (INT_MAX - 2) / 2)
                throw std::length_error("sequence exceeds the 32-bit alignment index range");
        omp_set_dynamic(0);
        omp_set_num_threads(o.threads);
        if (o.mode == "pair") run_pair(o, kernel, records);
        else run_matrix(o, kernel, records);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "editdist: " << e.what() << '\n';
        return 1;
    }
}
