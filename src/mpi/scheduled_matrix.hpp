#pragma once
#include "chps/dp.hpp"
#include "chps/fasta.hpp"
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <atomic>
#include <array>
#include <climits>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <thread>
#include <stdexcept>
#include <string>
#include <vector>

namespace chps::mpi {
struct Task {
    int i, j, slot;
    std::uint64_t order;
};

inline std::uint64_t hilbert_index(int side, int x, int y) {
    std::uint64_t index = 0;
    for (int scale = side / 2; scale; scale /= 2) {
        const int rx = (x & scale) != 0, ry = (y & scale) != 0;
        index += static_cast<std::uint64_t>(scale)*scale*((3*rx) ^ ry);
        if (!ry) {
            if (rx) {
                x = side - 1 - x;
                y = side - 1 - y;
            }
            std::swap(x,y);
        }
    }
    return index;
}

// Split the ordered tasks at the prefix nearest each target cumulative cost.
// Integer division avoids overflowing total*rank even when total fits uint64_t.
inline std::vector<long long> weighted_cuts(const std::vector<Task>& tasks,
                                             const std::vector<Record>& data, int ranks) {
    std::vector<std::uint64_t> prefix(tasks.size()+1,0);
    for (std::size_t i=0; i<tasks.size(); ++i) {
        const auto& t = tasks[i];
        const std::uint64_t a = data[t.i].bases.size(), b = data[t.j].bases.size();
        if (a && b>(std::numeric_limits<std::uint64_t>::max()-prefix[i])/a)
            throw std::length_error("weighted matrix exceeds 64-bit cell count");
        prefix[i+1] = prefix[i]+a*b;
    }
    std::vector<long long> cuts(ranks+1);
    cuts.back() = static_cast<long long>(tasks.size());
    const auto total = prefix.back();
    for (int rank=1; rank<ranks; ++rank) {
        if (!total) {
            cuts[rank]=cuts.back()*rank/ranks;
            continue;
        }
        const auto target = (total/ranks)*rank+(total%ranks)*rank/ranks;
        const auto upper = std::lower_bound(prefix.begin(),prefix.end(),target);
        auto offset = upper-prefix.begin();
        if (offset>0 && target-prefix[offset-1] < prefix[offset]-target) --offset;
        cuts[rank] = offset;
    }
    return cuts;
}

struct RankWork {
    unsigned long long pairs = 0, cells = 0, stolen = 0, cross_node = 0;
    double busy_seconds = 0;
};

struct Range {
    long long first = 0, last = 0;
    int owner = 0;
};

// A rank serves requests from its original contiguous block. Its workers take
// from the front, and remote thieves receive disjoint intervals from the back.
class VictimQueues {
    static constexpr int request_tag = 21900, reply_tag = 21901;
    int rank_, ranks_, grain_, gift_, victim_step_ = 1;
    Range own_, received_;
    bool waiting_ = false;
    std::vector<std::array<long long,2>> replies_;
    std::vector<MPI_Request> sends_;

    void progress() {
        MPI_Status status;
        int ready;
        for (;;) {
            MPI_Iprobe(MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,&ready,&status);
            if (!ready) break;
            const int source = status.MPI_SOURCE;
            MPI_Recv(nullptr,0,MPI_BYTE,source,request_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // A thief contacts this victim again only after receiving its reply.
            MPI_Wait(&sends_[source],MPI_STATUS_IGNORE);
            auto& reply = replies_[source];
            reply[1] = own_.last;
            own_.last = std::max(own_.first,own_.last-gift_);
            reply[0] = own_.last;
            MPI_Isend(reply.data(),2,MPI_LONG_LONG_INT,source,reply_tag,MPI_COMM_WORLD,&sends_[source]);
        }
        if (waiting_) {
            const int victim = (rank_+victim_step_)%ranks_;
            MPI_Iprobe(victim,reply_tag,MPI_COMM_WORLD,&ready,MPI_STATUS_IGNORE);
            if (ready) {
                std::array<long long,2> reply;
                MPI_Recv(reply.data(),2,MPI_LONG_LONG_INT,victim,reply_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                received_ = {reply[0],reply[1],victim};
                if (received_.first==received_.last) ++victim_step_;
                waiting_ = false;
            }
        }
    }

public:
    VictimQueues(long long count, int rank, int ranks, int grain, int team)
        : rank_(rank), ranks_(ranks), grain_(grain), gift_(grain * static_cast<long long>(team) > INT_MAX
            ? INT_MAX : grain*team), own_{count*rank/ranks,count*(rank+1)/ranks,rank},
          replies_(ranks), sends_(ranks,MPI_REQUEST_NULL) {}

    bool next(Range& chunk) {
        // All MPI calls from the OpenMP team pass this serialized region.
        bool keep_polling;
        #pragma omp critical(chps_mpi_steal)
        {
            progress();
            auto& source = own_.first<own_.last ? own_ : received_;
            chunk = {source.first,std::min(source.last,source.first+grain_),source.owner};
            source.first = chunk.last;
            keep_polling = chunk.first!=chunk.last || victim_step_<ranks_;
            if (chunk.first==chunk.last && victim_step_<ranks_ && !waiting_) {
                MPI_Request request;
                const int victim=(rank_+victim_step_)%ranks_;
                MPI_Isend(nullptr,0,MPI_BYTE,victim,request_tag,MPI_COMM_WORLD,&request);
                MPI_Request_free(&request);
                waiting_=true;
            }
        }
        return keep_polling;
    }

    void finish() {
        // A drained rank must continue to return empty replies until every rank
        // has finished claiming, computing and receiving its outstanding grants.
        MPI_Request barrier;
        MPI_Ibarrier(MPI_COMM_WORLD,&barrier);
        int done = 0;
        while (!done) {
            progress();
            MPI_Test(&barrier,&done,MPI_STATUS_IGNORE);
        }
        MPI_Waitall(ranks_,sends_.data(),MPI_STATUSES_IGNORE);
    }
};

inline std::vector<int> scheduled_matrix(const std::vector<Record>& data, const std::string& schedule,
                                          int rank, int ranks, int team, int grain,
                                          const std::vector<int>& host, RankWork& work) {
    const auto n = static_cast<long long>(data.size());
    const auto count = n * (n-1) / 2;
    if (count > INT_MAX) throw std::length_error("triangle exceeds MPI reduction count");
    const bool weighted = schedule=="hilbert-weighted";
    const bool curve = schedule=="hilbert" || weighted;
    const bool cyclic = schedule=="cyclic", steal = schedule=="steal";
    int side = 1;
    while (side < n) side *= 2;
    std::vector<Task> tasks;
    tasks.reserve(count);
    for (int i=0; i<n; ++i) for (int j=i+1; j<n; ++j) {
        const int p=static_cast<int>(tasks.size());
        tasks.push_back({i,j,p,curve ? hilbert_index(side,i,j) : 0});
    }
    if (curve) std::sort(tasks.begin(),tasks.end(),[](const Task& a,const Task& b){return a.order<b.order;});
    std::vector<int> scores(count,-1), result(rank ? 0 : count);
    std::vector<RankWork> per_worker(team);
    std::atomic<long long> next{0};
    std::exception_ptr failure;
    VictimQueues queues(count,rank,ranks,grain,team);
    const auto cuts = weighted ? weighted_cuts(tasks,data,ranks) : std::vector<long long>{};
    const auto begin = weighted ? cuts[rank] : cyclic ? rank : count*rank/ranks;
    const auto stride = cyclic ? ranks : 1;
    const auto size = weighted ? cuts[rank+1]-begin :
                      cyclic ? (count+ranks-1-rank)/ranks : count*(rank+1)/ranks-begin;
    #pragma omp parallel num_threads(team)
    {
        auto& local = per_worker[omp_get_thread_num()];
        for (;;) {
            Range chunk;
            if (steal) {
                if (!queues.next(chunk)) break;
                if (chunk.first==chunk.last) {
                    std::this_thread::yield();
                    continue;
                }
            } else {
                const auto pos = next.fetch_add(grain,std::memory_order_relaxed);
                if (pos>=size) break;
                chunk = {begin+pos*stride,begin+std::min(pos+grain,size)*stride,rank};
            }
            const double start = omp_get_wtime();
            for (auto i=chunk.first; i<chunk.last; i+=stride) {
                const auto& t=tasks[i];
                try {
                    scores[t.slot]=nvpd_omp(data[t.i].bases,data[t.j].bases).distance;
                    ++local.pairs;
                    local.cells+=static_cast<unsigned long long>(data[t.i].bases.size())*data[t.j].bases.size();
                    local.stolen+=chunk.owner!=rank;
                    local.cross_node+=host[chunk.owner]!=host[rank];
                } catch (...) {
                    #pragma omp critical(chps_schedule_error)
                    { if (!failure) failure=std::current_exception(); }
                }
            }
            local.busy_seconds+=omp_get_wtime()-start;
        }
    }
    if (steal) queues.finish();
    if (failure) std::rethrow_exception(failure);
    for (const auto& r : per_worker) {
        work.pairs+=r.pairs; work.cells+=r.cells; work.stolen+=r.stolen; work.cross_node+=r.cross_node;
        work.busy_seconds+=r.busy_seconds;
    }
    MPI_Reduce(scores.data(),result.data(),static_cast<int>(count),MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    if (!rank && std::find(result.begin(),result.end(),-1)!=result.end())
        throw std::runtime_error("unassigned matrix pair");
    return result;
}

inline std::vector<int> node_leaders(int rank, int ranks) {
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
    int leader=rank;
    MPI_Bcast(&leader,1,MPI_INT,0,node);
    MPI_Comm_free(&node);
    std::vector<int> all(ranks);
    MPI_Allgather(&leader,1,MPI_INT,all.data(),1,MPI_INT,MPI_COMM_WORLD);
    return all;
}

inline void rank_diagnostics(const RankWork& work, int rank, int ranks, int trial) {
    unsigned long long counters[]{work.pairs,work.cells,work.stolen,work.cross_node};
    std::vector<unsigned long long> all(rank ? 0 : 4*ranks);
    std::vector<double> times(rank ? 0 : ranks);
    MPI_Gather(counters,4,MPI_UNSIGNED_LONG_LONG,all.data(),4,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
    MPI_Gather(&work.busy_seconds,1,MPI_DOUBLE,times.data(),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (!rank) for (int i=0;i<ranks;++i)
        std::cerr << std::setprecision(12) << "{\"rank\":" << i << ",\"trial\":" << trial
                  << ",\"pairs\":" << all[4*i] << ",\"cells\":" << all[4*i+1]
                  << ",\"stolen_pairs\":" << all[4*i+2] << ",\"cross_node_pairs\":" << all[4*i+3]
                  << ",\"busy_seconds\":" << times[i] << "}\n";
}
} // namespace chps::mpi
