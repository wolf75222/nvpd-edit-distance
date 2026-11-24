# CPU profiling on ROMEO

The driver compares cost-only NvPD kernels on the same deterministic sequence pairs. It computes each expected distance with the two-row scalar recurrence before accepting a result. Every variant gets the same compiler, flags, input seed and OpenMP settings. Each timing trial starts a new process and warms the same kernel once. The clock surrounds the complete second call, including match-index construction, allocation, initialization and release. Process startup, input generation, the oracle and CPU-affinity reporting fall outside the interval.

Check out this repository with its history. The script retrieves the original CPU sources at `f08b4f6` and the vectorized version at `4724f68`. It loads the compact-index experiment from the tag `experiment/compact-index`. The script checks the source digests for all three variants, including explicitly staged sources. The experiment uses Linux, a C++20 compiler with OpenMP, Python 3, Slurm and a recent perf executable with FIFO control. Set `PERF_BIN` if perf is outside the job's PATH.

```bash
sbatch --account="$ROMEO_PROJECT" bench/profile_cpu.sbatch baseline
sbatch --account="$ROMEO_PROJECT" bench/profile_cpu.sbatch v1
sbatch --account="$ROMEO_PROJECT" bench/profile_cpu.sbatch v2
```

Each job writes to `results/profile/cpu/<job>/`. The batch script accepts optional `BASELINE_DIR`, `V1_SRC` and `V2_SRC` for clusters that stage source files without Git. `source-*.sha256` identifies the compiled translation units. The GCC vectorization report, flags, cache topology and Slurm allocation accompany the measurements.

## Designs

- The A/B experiments compare the initial branch-based kernel against branch-free SIMD (V1), then V1 against the 8-bit match-gap candidate (V2). Within each trial block the runner randomizes configuration order with a fixed seed. It alternates variant order between successive blocks. Three trials cover 1, 8 and 16 cores on each shape. For each pair of matching trial IDs, divide the initial-variant time by the comparison-variant time. Report the median of those three paired ratios.
- The strong-scaling case fixes the matrix at 20,000 × 20,000. The weak-scaling case fixes the row count at 4,096 and the column count at 4,096 times the worker count. Every worker receives exactly 16,777,216 interior cells. Both use 1, 2, 4, 8 and 16 workers. The hotter or colder per-thread cache working set can affect measured efficiency.
- The inputs contain uniform ACGT. Each overlapping position has a 5% substitution probability and the longer sequence has a random suffix. Cases are deterministic algorithm probes. The course corpus is generated separately by `scripts/synthetic.py`. In addition to the square shape, the 256 × 1,000,000 case probes long rows. The 64 × 4,194,304 case deliberately exceeds a 32 MiB EPYC CCD cache with its two integer rows and integer match index. Each matrix fits in RAM.

## Counters and memory

Perf counts user-space events for the parent and inherited OpenMP workers. The driver sends enable and disable through acknowledged FIFOs around a single warmed NvPD call. Measurement includes the clock calls and the small gate-handshake overhead. The events span three separate executions to limit multiplexing. Each pass repeats twice on 1 and 8 cores. `running_pct` comes from perf. The generic cache-reference and cache-miss counters retain their Linux event names. The L1 data and data-TLB events are separate counters. There is no DRAM byte measurement.

`dp_bytes` counts two integer score rows. `mit_bytes` counts the steady match-index bytes for the four letters. `workspace_peak_bytes` models the maximum simultaneously live score/index buffers. V2 first holds the integer index and its compact copy. It releases one before allocating score rows. This model omits the strings, small metadata, allocator and runtime. `peak_rss_kib` is the process high-water mark and includes the warm-up. `max_gap` verifies which V2 index path runs: all recorded cases use one-byte gaps. `cpu_ids` lists the logical CPUs observed by the OpenMP team after the timed call.
