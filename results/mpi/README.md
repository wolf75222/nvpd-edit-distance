# ROMEO MPI strong scaling

Job 838226 completed on `romeo-c009` and `romeo-c039`. Each step used 24 ranks and eight OpenMP workers per active node. The Linux placement checks counted 192 distinct physical cores on one node and 384 on two. Each cell count and input SHA-256 stayed fixed across the two layouts.

| Input | One node, s | Two nodes, s | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| One pair, 20,000 × 20,037 | 0.0911 | 0.9903 | 0.092× | 4.6% |
| One pair, 1,024 × 1,000,003 | 0.0148 | 0.0526 | 0.281× | 14.1% |
| 64 synthetic sequences, 2,016 pairs | 0.7643 | 0.5278 | 1.448× | 72.4% |
| 64 phiX174 genomes, 2,016 pairs | 0.1624 | 0.1146 | 1.417× | 70.8% |

Times are medians of three trials after one warmup. The single-pair algorithm exchanges two integers per rank on each row. On the shorter, wider pair it made 1,024 collective calls. On the nearly square pair it made 20,000. Adding the second node slowed both pairs. In matrix mode, each worker evaluated independent pairs before the ranks reduced the result. The two matrix benchmarks completed faster on two nodes.

[Trial data](838226/trials.csv) include the timed region, assigned workers, sampled physical cores, input digests and edlib checks. The [summary](838226/summary.csv) uses only job 838226. The [MPI smoke tests](838226/mpi-tests.txt) cover empty inputs, uneven stripes, distant matches, excess ranks and IUPAC symbols.

Job 838221 stopped at MPI startup. Job 838222 completed the two single-pair comparisons and the separate Kokkos CPU check. Its 192-rank-per-node matrix launch exceeded the 24 MPI slots per node provided by that allocation. Its pair-only trials are kept separately from the completed campaign. The [Slurm accounting](jobs.tsv) records all three states.

The [source digests](838226/source.sha256) identify the files used in the final job. The published batch script subsequently restricted Slurm environment logging to the allocation fields displayed in the results. See the [MPI protocol](../../bench/mpi.md) for commands and timing boundaries.

The [weak-scaling study](weak/README.md) doubles the actual DP work at constant cells per physical core and reports throughput and weak efficiency separately.
