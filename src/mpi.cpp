#include "chps/dp.hpp"
#include "chps/fasta.hpp"
#include "dp_internal.hpp"
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <charconv>
#include <climits>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#ifdef __linux__
#include <sched.h>
#endif

namespace {
struct Options {
    std::string mode, input;
    int threads = 1, trials = 1, warmups = 0;
};

Options parse(int argc, char** argv) {
    if (argc < 2) throw std::invalid_argument("expected pair or matrix");
    Options o{argv[1], {}};
    if (o.mode != "pair" && o.mode != "matrix") throw std::invalid_argument("unknown mode");
    for (int i = 2; i < argc; i += 2) {
        if (i + 1 >= argc) throw std::invalid_argument("missing argument value");
        std::string_view key = argv[i], value = argv[i+1];
        if (key == "--input") { o.input = value; continue; }
        int n = 0;
        auto [end, ec] = std::from_chars(value.data(), value.data() + value.size(), n);
        if (ec != std::errc{} || end != value.data() + value.size())
            throw std::invalid_argument("expected integer for " + std::string(key));
        if (key == "--threads" && n > 0) o.threads = n;
        else if (key == "--trials" && n > 0) o.trials = n;
        else if (key == "--warmups" && n >= 0) o.warmups = n;
        else throw std::invalid_argument("unknown option or invalid value: " + std::string(key));
    }
    if (o.input.empty()) throw std::invalid_argument("--input is required");
    return o;
}

void broadcast_string(std::string& text, int rank) {
    if (text.size() > INT_MAX) throw std::length_error("MPI message exceeds INT_MAX bytes");
    int length = static_cast<int>(text.size());
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank) text.resize(length);
    MPI_Bcast(text.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
}

std::vector<chps::Record> input(const Options& o, int rank) {
    std::vector<chps::Record> records;
    if (!rank) records = chps::read_fasta(o.input);
    if (records.size() > INT_MAX) throw std::length_error("too many FASTA records");
    int count = static_cast<int>(records.size());
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (o.mode == "pair" && count != 2) throw std::invalid_argument("pair requires two records");
    if (rank) records.resize(count);
    for (auto& r : records) {
        broadcast_string(r.id, rank);
        broadcast_string(r.bases, rank);
    }
    return records;
}

// All coordinates are one-based in B. Rank r owns columns starts[r]+1..starts[r+1].
int distributed_pair(std::string_view a, std::string_view b, int rank, int ranks, int team,
                     int& workers) {
    chps::detail::validate(a, b, false);
    if (a.size() > b.size()) std::swap(a,b);
    const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    if (!m) { workers = 0; return n; }
    std::vector<int> starts(ranks + 1);
    for (int r = 0; r <= ranks; ++r)
        starts[r] = (n / ranks) * r + std::min(n % ranks, r);
    const int start = starts[rank], width = starts[rank+1] - start;
    workers = std::min(team, width);
    const auto local = b.substr(start, width);
    const chps::detail::MatchIndex matches(a, local);
    std::array<int, 256> preceding{}, slot{};
    for (int j = 0; j < start; ++j) preceding[static_cast<unsigned char>(b[j])] = j + 1;
    // Only one match before the stripe can contribute to each symbol. Its
    // predecessor is the owning rank's last-match value, or its left neighbor's
    // right edge when the match lies in the first column of its stripe.
    for (std::size_t s = 0; s < slot.size(); ++s) {
        const int k = preceding[s];
        if (!k) continue;
        const int owner = static_cast<int>(std::upper_bound(starts.begin(), starts.end(), k-1)
                                           - starts.begin()) - 1;
        slot[s] = k == starts[owner] + 1 ? 2 * (owner - 1) : 2 * owner + 1;
    }
    std::vector<int> edges(2 * static_cast<std::size_t>(ranks)), prev(width+1), curr(width+1);
    std::iota(prev.begin(), prev.end(), start);
    for (int i = 1; i <= m; ++i) {
        const auto c = static_cast<unsigned char>(a[i-1]);
        const int* last = matches.entries.data() + matches.symbol[c] * (local.size()+1);
        const int k_end = last[width];
        const int send[] = {prev[width], k_end > 1 ? prev[k_end-1] : 0};
        MPI_Allgather(send, 2, MPI_INT, edges.data(), 2, MPI_INT, MPI_COMM_WORLD);
        prev[0] = rank ? edges[2*(rank-1)] : i-1;
        const int outer = slot[c] >= 0 ? edges[slot[c]] : i-1;
        const int k_left = preceding[c];
        #pragma omp parallel for simd schedule(static) if(team > 1)
        for (int x = 1; x <= width; ++x) {
            const int k = last[x];
            const int upper = 1 + std::min(prev[x-1], prev[x]);
            const int via = k ? prev[k-1] + (x-k) :
                            k_left ? outer + (start+x-k_left) : upper;
            curr[x] = std::min(upper, via);
        }
        std::swap(prev,curr);
    }
    int result = prev[width];
    MPI_Bcast(&result, 1, MPI_INT, std::min(ranks,n)-1, MPI_COMM_WORLD);
    return result;
}

std::vector<int> batch(const std::vector<chps::Record>& data, int rank, int ranks,
                       int team, int& workers) {
    const std::size_t count = data.size() * (data.size()-1) / 2;
    if (count > INT_MAX) throw std::length_error("cohort exceeds MPI reduction count limit");
    std::vector<int> local(count), result(rank ? 0 : count);
    std::vector<std::pair<std::size_t,std::size_t>> jobs;
    std::size_t p = 0;
    for (std::size_t i = 0; i < data.size(); ++i) {
        for (std::size_t j = i+1; j < data.size(); ++j, ++p) {
            if (p % ranks == static_cast<std::size_t>(rank)) jobs.emplace_back(i,j);
        }
    }
    workers = static_cast<int>(std::min<std::size_t>(team,jobs.size()));
    // Each worker owns a whole pair. Nested OpenMP is disabled, so the inner
    // NvPD uses one worker and its SIMD columns without oversubscribing the node.
    std::exception_ptr failure;
    #pragma omp parallel for schedule(dynamic,1)
    for (std::size_t task = 0; task < jobs.size(); ++task) {
        const auto [i,j] = jobs[task];
        try {
            local[rank + task*ranks] = chps::nvpd_omp(data[i].bases,data[j].bases).distance;
        } catch (...) {
            #pragma omp critical
            { if (!failure) failure = std::current_exception(); }
        }
    }
    if (failure) std::rethrow_exception(failure);
    MPI_Reduce(local.data(), result.data(), static_cast<int>(count), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return result;
}

void matrix_output(const std::vector<chps::Record>& data, const std::vector<int>& triangle) {
    std::cout << "id";
    for (const auto& r : data) std::cout << '\t' << r.id;
    std::cout << '\n';
    const auto n = data.size();
    for (std::size_t i = 0; i < n; ++i) {
        std::cout << data[i].id;
        for (std::size_t j = 0; j < n; ++j) {
            const auto [lo,hi] = std::minmax(i,j);
            std::cout << '\t' << (i == j ? 0 : triangle[lo*(2*n-lo-1)/2 + hi-lo-1]);
        }
        std::cout << '\n';
    }
}

// Sample each OpenMP worker; collapse SMT siblings by package/core on each host.
int placement(int threads, int& team) {
    std::vector<int> cpus(threads,-1);
    #pragma omp parallel
    {
        #pragma omp single
        team = omp_get_num_threads();
#ifdef __linux__
        cpus[omp_get_thread_num()] = sched_getcpu();
#endif
    }
#ifdef __linux__
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&node);
    int rank, size;
    MPI_Comm_rank(node,&rank);
    MPI_Comm_size(node,&size);
    std::vector<int> ids(rank ? 0 : threads*size);
    MPI_Gather(cpus.data(),threads,MPI_INT,ids.data(),threads,MPI_INT,0,node);
    MPI_Comm_free(&node);
    std::set<std::pair<int,int>> physical;
    if (!rank) for (int cpu : ids) {
        if (cpu < 0) continue;
        const auto dir = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/";
        int socket, core;
        if (!(std::ifstream(dir+"physical_package_id") >> socket) || !(std::ifstream(dir+"core_id") >> core))
            throw std::runtime_error("cannot resolve the physical core of an MPI worker");
        physical.emplace(socket,core);
    }
    const int local = static_cast<int>(physical.size());
    int global;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    return global;
#else
    return -1;
#endif
}
} // namespace

int main(int argc, char** argv) {
    int support, rank, ranks;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &support);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try {
        if (support < MPI_THREAD_FUNNELED) throw std::runtime_error("MPI_THREAD_FUNNELED required");
        const auto o = parse(argc,argv);
        const auto records = input(o,rank);
        omp_set_dynamic(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(o.threads);
        int team = 1;
        const int physical = placement(o.threads,team);
        int distance = 0;
        std::vector<int> matrix;
        int first_distance = 0;
        std::vector<int> first_matrix;
        for (int trial = -o.warmups; trial < o.trials; ++trial) {
            MPI_Barrier(MPI_COMM_WORLD);
            const double start = MPI_Wtime();
            int workers = 0;
            if (o.mode == "pair")
                distance = distributed_pair(records[0].bases, records[1].bases, rank, ranks, team, workers);
            else matrix = batch(records, rank, ranks, team, workers);
            double elapsed = MPI_Wtime() - start, slowest = 0;
            MPI_Reduce(&elapsed, &slowest, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            int total_workers = 0;
            MPI_Reduce(&workers, &total_workers, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            if (!rank) {
                if (trial == -o.warmups) { first_distance = distance; first_matrix = matrix; }
                else if (first_distance != distance || first_matrix != matrix)
                    throw std::runtime_error("repeated MPI distances disagree");
            }
            if (!rank && trial >= 0)
                std::cerr << std::setprecision(12) << "{\"mode\":\"" << o.mode << "\",\"trial\":" << trial+1
                          << ",\"ranks\":" << ranks << ",\"threads_per_rank\":" << team
                          << ",\"sampled_physical_cores\":" << physical
                          << ",\"assigned_workers\":" << total_workers << ",\"seconds\":" << slowest << "}\n";
        }
        if (!rank) {
            if (o.mode == "matrix") matrix_output(records,matrix);
            else std::cout << "{\"algorithm\":\"nvpd-mpi\",\"distance\":" << distance << "}\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "editdist-mpi rank " << rank << ": " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
}
