# Distributed NvPD

Build with an MPI implementation that supports `MPI_THREAD_FUNNELED`, a C++20 compiler and OpenMP.

```bash
cmake -S . -B build-mpi -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON
cmake --build build-mpi --target editdist-mpi editdist -j
mpiexec -n 4 build-mpi/editdist-mpi pair --input pair.fa --threads 2
mpiexec -n 4 build-mpi/editdist-mpi matrix --input cohort.fa --threads 2 > matrix.tsv
python3 scripts/tree.py matrix.tsv --method nj --out tree.nwk
```

Rank zero reads the FASTA file and broadcasts its records. The pair command writes a JSON distance. The matrix command writes a labeled square TSV. Timing records go to stderr. The options `--warmups N` and `--trials N` repeat the calculation on the same input. Every repetition must return identical results.

## One pair

NvPD advances one row of the dynamic program at a time. The shorter sequence supplies the rows. MPI ranks divide the columns into contiguous stripes, including empty stripes if the input has fewer columns than ranks. Each rank stores a local match index and two score rows. OpenMP threads compute the columns of that rank.

For row character *c*, the previous occurrence in a target prefix can lie several ranks to the left. The target stays fixed throughout the alignment, so each rank precomputes where that predecessor cost will come from for each symbol. At each row, an MPI allgather exchanges two integers per rank. The first is the previous row's right-edge cost. The second is the predecessor cost for the last local match of *c*. A match at the stripe's first column takes its predecessor from the preceding rank's right edge. Column zero uses the row number on the global boundary. If *c* has no previous match, the ordinary deletion and substitution candidates suffice.

This schedule visits all *mn* cells and performs one collective per row. A rank owns at most `ceil(n / ranks)` score columns. All ranks hold a copy of both input sequences. Adding nodes reduces the number of local columns and introduces inter-node communication on every row.

## Many pairs

For a cohort, the ranks distribute the upper triangular pairs in round-robin order. Within a rank each OpenMP worker takes a complete pair and evaluates the SIMD NvPD recurrence. A single MPI reduction assembles the exact distances on rank zero. The output diagonal is zero. This mode measures throughput across independent pairs.

## ROMEO

```bash
sbatch --account="$ROMEO_PROJECT" --partition=instant bench/mpi.sbatch
```

The job requests two exclusive x64cpu nodes with 192 physical cores and 16 GiB per node. Both modes use 24 ranks per active node and eight OpenMP workers per rank. The one-node and two-node measurements run in the same reservation. Each input gets one warmup and three timed calls in each layout. The benchmark excludes file input, initial input broadcast, MPI startup and final text output. Each trial records the maximum elapsed time across ranks, including local setup, computation and MPI transfers within the backend. On Linux, the executable samples each OpenMP worker's CPU and identifies its physical core from sysfs before the timed region.

The job checks MPI on two hosts before running the study. It uses the Slurm PMIx launcher if the probe succeeds, or the allocated hosts through Open MPI if needed. The optional `KOKKOS_SOURCE` environment variable requests a separate Kokkos CPU check before the MPI phase.

The measured pairs have dimensions 20,000 × 20,037 and 1,024 × 1,000,003. The synthetic cohort contains 64 sequences of 12,000 bases. If `data/phix174/cohort.fa` is present, the job also evaluates all 2,016 phage pairs. Generate that FASTA with the project's phage dataset script. Each final result is compared with edlib. For a matrix this means comparing the SHA-256 of the complete labeled matrix. The smoke tests also compare small pairs to classical sequential DP.

Results appear under `results/mpi/<job>/`. For each input, the strong speedup is the one-node median divided by the two-node median. Efficiency divides that ratio by two. The column `assigned_workers` counts workers assigned DP columns or cohort tasks. The column `sampled_physical_cores` counts distinct observed physical-core placements. Both are distinct from a time-averaged CPU utilization measurement.
