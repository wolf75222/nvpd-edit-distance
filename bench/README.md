# Benchmarks on ROMEO

Submit from the project root. Both jobs build in a temporary directory on the compute node. The CPU job requests eight cores. The GPU job requests four host cores and one Hopper GPU. Both request 8 GiB for 30 minutes. The build needs the two dependency archives in `deps/` or network access to the sources configured by CMake.

```bash
python3 scripts/datasets/synthetic.py --out data/synthetic
sbatch --account="$ROMEO_PROJECT" bench/campaign/cpu.sbatch
sbatch --account="$ROMEO_PROJECT" bench/campaign/gpu.sbatch
python3 bench/campaign/summarize.py results/cpu-*.csv results/gpu-*.csv > results/summary.csv
```

The default run uses the first pair at each error rate for 5,000 and 20,000 nt. Each configuration runs three times in a fresh process. OpenMP runs on 1, 2, 4 and 8 cores with fixed thread placement. The Slurm log records the CPU, CUDA device and compiler configuration. A CSV row records the job and host, plus SHA-256 digests of the executable and input FASTA. The archived campaign includes the [submission commands](../results/romeo/submissions.tsv) and [kernel source checksums](../results/romeo/sources.sha256).

Use `--profile full` to select every pair. Split long jobs by length, rate, replicate or backend. Adjust the Slurm wall time to the measured runtime. At 1,000,000 nt, the sequential DP reference runs once per rate on replicate 1. Pair-specific speedups stay empty for replicates without a sequential measurement. Example of a single shard:

```bash
sbatch --account="$ROMEO_PROJECT" --time=02:00:00 bench/campaign/cpu.sbatch \
  --profile full --lengths 100000 --rates 0.02 --replicates 1 2 \
  --algorithms seq nvpd-omp --threads 1 4 8 --trials 1
```

For offline compute nodes, stage the source directories, the two archives in `deps/`, and `data/synthetic/` on ROMEO before submitting. Stage the cohort FASTA as well when requesting a matrix.

The Python runner also accepts a prebuilt executable. Use `python3 bench/campaign/run.py --help` for the arguments. It flushes each trial to disk. A failed or timed-out subprocess leaves a diagnostic row. The job exits with an error if any trial fails.

## Reading the results

The summary compares the exact distances for every matching FASTA digest and stops on disagreement. For a fixed campaign, pass `--plan submissions.tsv` with columns `job` and `expected_trials` to detect missing rows or jobs before aggregation. Inspect the Slurm exit state for timed-out allocations. The summary reports the median, minimum and maximum of alignment times. `median_wall_seconds` includes process launch, FASTA input, alignment and JSON output. `max_peak_rss_kib` covers the whole host process. `max_device_bytes` counts explicit device buffer allocations.

For lengths *m* and *n* and median alignment time *T*, `equivalent_gcups` equals *mn* / (10⁹ *T*). This quantity measures equivalent dense-DP throughput for all algorithms. `dense_gcups` applies to the four implementations that visit all interior cells. The library backends use different work units.

`speedup_vs_cpu_seq` divides the sequential CPU time by the backend time for the same pair. For a CPU job, the summary prefers its own reference. Otherwise it uses the median of the supplied CPU sequential trials for that pair. The `sequential_reference_jobs` field identifies the contributing job IDs. The GPU speedup compares two machines. For OpenMP, `omp_efficiency_vs_cpu_seq` divides this speedup by the number of threads. `omp_scaling` compares a backend to its own single-thread trial on the same host and job. `omp_efficiency` divides this ratio by the thread count.

## Cohort matrix and trees

Both Slurm jobs run the native differential tests before timing the synthetic pairs. The GPU job exercises both CUDA storage modes and replays the traceback.

Set `COHORT_FASTA` to run both library backends over every unordered pair. Set `COHORT_VERIFY_WITH=nvpd-omp` or, on a GPU job, `COHORT_VERIFY_WITH=nvpd-cuda`, to check the same matrix with NvPD. Supply an optional TSV with `id` and `group` columns via `COHORT_GROUPS`. The workflow compares both matrices byte for byte and records each process wall time. It constructs UPGMA and neighbor-joining trees from the verified matrix. Files go to `results/cohort-<job>` independently of the synthetic trials.

```bash
COHORT_FASTA=data/cohort.fa COHORT_GROUPS=data/groups.tsv \
    sbatch --account="$ROMEO_PROJECT" bench/campaign/cpu.sbatch
```

## Profilage et architectures

| Expérience | Protocole | Mesures |
|---|---|---|
| SIMD, caches, RSS, échelle forte et faible | [OpenMP](cpu/openmp/profile_cpu.md) | [Compteurs et temps](../results/profile/cpu/) |
| CUDA, tailles de tuiles, Nsight Systems | [Script Slurm](gpu/cuda/profile_gpu.sbatch) | [Appels, noyaux et gains](../results/profile/gpu/) |
| Deux nœuds complets | [MPI et OpenMP](mpi/mpi.md) | [Résultats à 192 et 384 cœurs](../results/mpi/) |
| Allocations GPU par lot | [Protocole Slurm](gpu/allocation/allocator.sbatch) | [Réutilisation des tampons](../results/allocator/) |
| Portabilité CPU/GPU | [Kokkos](kokkos/kokkos.md) | [Comparaison aux versions natives](../results/kokkos/) |
| Chromosome 21 | [Script Slurm](human/human.sbatch) | [Fenêtres et contrôle entier](../results/human/) |
| Méthodes en bande et bibliothèques | [Comparaison locale](cpu/alternatives/alternatives.py) | [Temps et empreintes](../results/alternatives/) |
| Fronts sur GH200 | [eWFA-GPU](gpu/ewfa/ewfa.sbatch) | [Validation des scores et des traces](../results/ewfa/) |
| Localité sur GH200 | [Secteurs par warp](gpu/locality/gpu_sectors.cu) | [Demandes de secteurs par case](../results/locality/) |
| Passage à l'échelle faible sur deux nœuds | [Paires et lots MPI](mpi/mpi_weak.py) | [Essais à charge par cœur constante](../results/mpi/weak/) |
| Quinze paires à un million | [Deux versions OpenMP](cpu/openmp/million/complete_1m.sbatch) | [Temps des quinze entrées](../results/romeo/completion-1m/) |
| Chromosomes entiers | [Limites de temps et d'espace](human/whole_chromosomes.py) | [Mesures et bornes](../results/human/whole/) |
| Trace sur processeur | [Matrice et Hirschberg](cpu/traceback/traceback.py) | [Temps et RSS](../results/traceback/) |
