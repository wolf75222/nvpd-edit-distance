#include "chps/dp.hpp"
#include "chps/fasta.hpp"
#include "dp_internal.hpp"
#ifdef CHPS_KOKKOS
#include "chps/kokkos.hpp"
#include <optional>
#endif
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
#include <sys/resource.h>
#ifdef __linux__
#include <sched.h>
#endif

#include "scheduled_matrix.hpp"

namespace {
struct Options {
    std::string mode, input;
    int threads = 1, trials = 1, warmups = 0, pipeline_rows = 0;
    bool report_rss = false;
    std::string schedule;
    bool kokkos_workers = false;
    int grain = 1;
};

Options parse(int argc, char** argv) {
    if (argc < 2) throw std::invalid_argument("expected pair, pairs or matrix");
    Options o{argv[1], {}};
    if (o.mode != "pair" && o.mode != "pairs" && o.mode != "matrix") throw std::invalid_argument("unknown mode");
    for (int i = 2; i < argc; i += 2) {
        if (i + 1 >= argc) throw std::invalid_argument("missing argument value");
        std::string_view key = argv[i], value = argv[i+1];
        if (key == "--input") { o.input = value; continue; }
        if (key == "--workers" && o.mode == "matrix") {
            if (value != "openmp" && value != "kokkos") throw std::invalid_argument("expected openmp or kokkos");
            o.kokkos_workers = value == "kokkos";
            continue;
        }
        if (key == "--schedule" && o.mode == "matrix") {
            if (value != "static" && value != "cyclic" && value != "hilbert" && value != "hilbert-weighted" && value != "steal")
                throw std::invalid_argument("expected static, cyclic, hilbert, hilbert-weighted or steal");
            o.schedule = value;
            continue;
        }
        int n = 0;
        auto [end, ec] = std::from_chars(value.data(), value.data() + value.size(), n);
        if (ec != std::errc{} || end != value.data() + value.size())
            throw std::invalid_argument("expected integer for " + std::string(key));
        if (key == "--threads" && n > 0) o.threads = n;
        else if (key == "--grain" && o.mode == "matrix" && n > 0) o.grain = n;
        else if (key == "--pipeline-rows" && o.mode == "pair" && n > 0) o.pipeline_rows = n;
        else if (key == "--report-rss" && n == 1) o.report_rss = true;
        else if (key == "--trials" && n > 0) o.trials = n;
        else if (key == "--warmups" && n >= 0) o.warmups = n;
        else throw std::invalid_argument("unknown option or invalid value: " + std::string(key));
    }
    if (o.kokkos_workers && o.schedule.empty()) throw std::invalid_argument("Kokkos workers require --schedule");
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
    if (o.mode == "pairs" && count % 2) throw std::invalid_argument("pairs requires an even record count");
    if (rank) records.resize(count);
    for (auto& r : records) {
        broadcast_string(r.id, rank);
        broadcast_string(r.bases, rank);
    }
    return records;
}

// Allgather recurrence, with the global last match before this stripe.
int allgather_stripe(std::string_view a, std::string_view b, const std::vector<int>& starts,
                     const chps::detail::MatchIndex& matches, int rank, int team) {
    const int m = static_cast<int>(a.size()), ranks = static_cast<int>(starts.size() - 1);
    const int start = starts[rank], width = starts[rank+1] - start;
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
        const int* last = matches.entries.data() + matches.symbol[c] * (static_cast<std::size_t>(width)+1);
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
    return prev[width];
}

// Each message contains the consecutive global left boundary scores of a
// rectangle. Its bottom boundary stays in this rank for the following block.
int pipeline_stripe(std::string_view a, int start, int width, int rank, int active_ranks,
                    int workers, int height, const chps::detail::MatchIndex& matches) {
    if (!width) return 0;
    const int m = static_cast<int>(a.size());
    height = std::min(height, m);
    const std::size_t stride = static_cast<std::size_t>(width) + 1;
    std::vector<int> even(stride), odd(stride), left(height), right(height);
    std::iota(even.begin(), even.end(), start);
    for (int done = 0; done < m;) {
        const int count = std::min(height, m - done);
        if (rank) MPI_Recv(left.data(), count, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else std::iota(left.begin(), left.begin() + count, done + 1);
        #pragma omp parallel num_threads(workers) if(workers > 1)
        {
            for (int t = 0; t < count; ++t) {
                const int i = done + t + 1;
                const int* prev = (i % 2 ? even : odd).data();
                int* curr = (i % 2 ? odd : even).data();
                const auto symbol = matches.symbol[static_cast<unsigned char>(a[i-1])];
                const int* last = matches.entries.data() + symbol * stride;
                const int edge = left[t];
                #pragma omp master
                curr[0] = edge;
                #pragma omp for simd schedule(static)
                for (int x = 1; x <= width; ++x) {
                    const int edit = 1 + std::min(prev[x], prev[x-1]);
                    const int from_edge = edge + x;
                    const int k = last[x];
                    const int match = k ? prev[k-1] + (x-k) : from_edge;
                    curr[x] = std::min({edit, from_edge, match});
                } // All readers finish before the older buffer is overwritten.
                #pragma omp master
                right[t] = curr[width];
            }
        }
        if (rank + 1 < active_ranks)
            MPI_Send(right.data(), count, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        done += count;
    }
    return (m % 2 ? odd : even)[width];
}

int distributed_pair(std::string_view a, std::string_view b, int rank, int ranks, int team,
                     int pipeline_rows, int& workers) {
    chps::detail::validate(a, b, false);
    if (a.size() > b.size()) std::swap(a,b);
    const int n = static_cast<int>(b.size());
    if (a.empty()) { workers = 0; return n; }
    std::vector<int> starts(static_cast<std::size_t>(ranks) + 1);
    for (int r = 0; r < ranks; ++r)
        starts[r] = (n / ranks) * r + std::min(n % ranks, r);
    starts[ranks] = n;
    const int start = starts[rank], width = starts[rank+1] - start;
    workers = std::min(team, width);
    const chps::detail::MatchIndex matches(a, b.substr(start, width));
    const int active_ranks = std::min(n, ranks);
    int result = pipeline_rows ? pipeline_stripe(a, start, width, rank, active_ranks, workers,
                                                 pipeline_rows, matches)
                               : allgather_stripe(a, b, starts, matches, rank, team);
    MPI_Bcast(&result, 1, MPI_INT, active_ranks - 1, MPI_COMM_WORLD);
    return result;
}

std::vector<int> batch(const std::vector<chps::Record>& data, int rank, int ranks,
                       int team, bool adjacent, int& workers) {
    const std::size_t count = adjacent ? data.size()/2 : data.size() * (data.size()-1) / 2;
    if (count > INT_MAX) throw std::length_error("cohort exceeds MPI reduction count limit");
    std::vector<int> local(count), result(rank ? 0 : count);
    std::vector<std::pair<std::size_t,std::size_t>> jobs;
    if (adjacent) {
        for (std::size_t p = rank; p < count; p += ranks) jobs.emplace_back(2*p, 2*p+1);
    } else {
        std::size_t p = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            for (std::size_t j = i+1; j < data.size(); ++j, ++p) {
                if (p % ranks == static_cast<std::size_t>(rank)) jobs.emplace_back(i,j);
            }
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
    int required = MPI_THREAD_FUNNELED;
    for (int i = 1; i+1 < argc; ++i)
        if (std::string_view(argv[i]) == "--schedule" && std::string_view(argv[i+1]) == "steal")
            required = MPI_THREAD_SERIALIZED;
    MPI_Init_thread(&argc, &argv, required, &support);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try {
        if (support < required) throw std::runtime_error("MPI thread support insufficient for chosen schedule");
        const auto o = parse(argc,argv);
        const auto records = input(o,rank);
        omp_set_dynamic(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(o.threads);
        int team = 1;
        const int physical = placement(o.threads,team);
#if defined(CHPS_KOKKOS) && defined(KOKKOS_ENABLE_OPENMP)
        std::optional<chps::KokkosSession> kokkos;
        if (o.kokkos_workers) kokkos.emplace(team);
#else
        if (o.kokkos_workers) throw std::invalid_argument("Kokkos workers require a Kokkos OpenMP build");
#endif
        const auto hosts = o.schedule.empty() ? std::vector<int>{} : chps::mpi::node_leaders(rank,ranks);
        int distance = 0;
        std::vector<int> matrix;
        int first_distance = 0;
        std::vector<int> first_matrix;
        for (int trial = -o.warmups; trial < o.trials; ++trial) {
            MPI_Barrier(MPI_COMM_WORLD);
            const double start = MPI_Wtime();
            int workers = 0;
            chps::mpi::RankWork work;
            if (o.mode == "pair")
                distance = distributed_pair(records[0].bases, records[1].bases, rank, ranks, team, o.pipeline_rows, workers);
            else if (!o.schedule.empty()) {
                matrix = chps::mpi::scheduled_matrix(records,o.schedule,rank,ranks,team,o.grain,hosts,work,o.kokkos_workers);
                workers = static_cast<int>(std::min<unsigned long long>(team,work.pairs));
            } else matrix = batch(records, rank, ranks, team, o.mode == "pairs", workers);
            double elapsed = MPI_Wtime() - start, slowest = 0;
            MPI_Reduce(&elapsed, &slowest, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            int total_workers = 0;
            MPI_Reduce(&workers, &total_workers, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            if (!rank) {
                if (trial == -o.warmups) { first_distance = distance; first_matrix = matrix; }
                else if (first_distance != distance || first_matrix != matrix)
                    throw std::runtime_error("repeated MPI distances disagree");
            }
            long max_rss_kib = 0;
            if (o.report_rss && trial >= 0) {
                rusage usage{};
                if (getrusage(RUSAGE_SELF,&usage)) throw std::runtime_error("getrusage failed");
                long rss_kib = usage.ru_maxrss;
#ifdef __APPLE__
                rss_kib /= 1024;
#endif
                MPI_Reduce(&rss_kib,&max_rss_kib,1,MPI_LONG,MPI_MAX,0,MPI_COMM_WORLD);
            }
            if (!o.schedule.empty() && trial >= 0)
                chps::mpi::rank_diagnostics(work,rank,ranks,trial+1);
            if (!rank && trial >= 0) {
                std::cerr << std::setprecision(12) << "{\"mode\":\"" << o.mode << "\",\"trial\":" << trial+1
                          << ",\"ranks\":" << ranks << ",\"threads_per_rank\":" << team
                          << ",\"sampled_physical_cores\":" << physical
                          << ",\"pipeline_rows\":" << o.pipeline_rows
                          << ",\"assigned_workers\":" << total_workers << ",\"seconds\":" << slowest;
                if (!o.schedule.empty()) std::cerr << ",\"worker_backend\":\"" << (o.kokkos_workers ? "kokkos" : "openmp") << '"';
                if (o.report_rss) std::cerr << ",\"max_rank_rss_kib\":" << max_rss_kib;
                std::cerr << "}\n";
            }
        }
        if (!rank) {
            if (o.mode == "matrix") matrix_output(records,matrix);
            else if (o.mode == "pairs") {
                std::cout << "id_a\tid_b\tdistance\n";
                for (std::size_t i = 0; i < matrix.size(); ++i)
                    std::cout << records[2*i].id << '\t' << records[2*i+1].id << '\t' << matrix[i] << '\n';
            }
            else std::cout << "{\"algorithm\":\"" << (o.pipeline_rows ? "nvpd-mpi-pipeline" : "nvpd-mpi")
                           << "\",\"distance\":" << distance << "}\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "editdist-mpi rank " << rank << ": " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
}
