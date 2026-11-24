# Full-node MPI weak scaling

Job 838444 reserved two exclusive 192-core EPYC nodes. Within the reservation, 24 MPI ranks per active node each used eight OpenMP workers. The placement check found 192 unique physical cores in every one-node trial and 384 in every two-node trial. The rank count rose from 24 to 48. Each case used one warmup followed by five timed calls per layout. Elapsed time is the maximum over ranks, including local setup, DP and inter-rank communication. The input read, broadcast, startup and output are outside the timed calls. [Placement and allocation](838444/environment.txt) · [Raw trials and input hashes](838444/trials.csv)

| Input, one → two nodes | DP cells per core | One node (s) | Two nodes (s) | Weak efficiency | Throughput ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1,024 rows, 49,152 → 98,304 columns | 262,144 | 0.005055 | 0.051055 | 9.90% | 0.198× |
| 1,024 rows, 786,432 → 1,572,864 columns | 4,194,304 | 0.013937 | 0.051977 | 26.81% | 0.536× |
| 85 → 120 sequences of 8,192 nt | 1,247,805,440 | 0.686621 | 0.893083 | 76.88% | 1.538× |

Times are medians. Weak efficiency is the one-node time divided by the two-node time. The throughput ratio compares total DP cells per second and equals twice the weak efficiency here. Each single-pair worker has exactly 256 or 4,096 columns and the same 1,024 rows on both layouts. The longer target extends the one-node target without changing its prefix. All 1,024 rows require a collective in each layout. The second node increases pair latency even with constant local DP work. [Summary](838444/summary.csv)

The larger cohort contains the initial 85 records. The number of unordered pairs grows from 3,570 to 7,140. Equal sequence lengths give identical mean cell counts per worker on both layouts. Matrix work is dynamically scheduled within each rank, so individual workers can complete different numbers of pairs. The full matrix throughput increases by 53.8%. Every pair score and every entry of each complete matrix agreed with edlib outside timing. The job also ran native and MPI differential tests. [Native tests](838444/native-tests.txt) · [MPI tests](838444/mpi-tests.txt)

Run from the repository root on ROMEO with `sbatch --account="$ROMEO_PROJECT" --partition=instant bench/mpi/mpi.sbatch weak`. Generation uses a fixed seed and writes temporary FASTA files. The script alternates which layout runs first across cases. Results are under `results/mpi/weak/$SLURM_JOB_ID/`. The intervals in the trial data are repeatability observations in one node allocation. [Generator and collector](../../../bench/mpi/mpi_weak.py) · [Measured source hashes](838444/source.sha256) · [Slurm accounting](838444/accounting.tsv)
