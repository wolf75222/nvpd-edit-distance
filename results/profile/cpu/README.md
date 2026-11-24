# CPU measurements on ROMEO

Run the three Slurm commands in [the CPU profiling protocol](../../../bench/profile_cpu.md) from the repository root to reproduce these experiments. The compiler was GCC 11.5.0 with `-std=c++20 -O3 -g -march=native -fopenmp -fno-omit-frame-pointer`. The jobs placed OpenMP workers on physical cores with close affinity, active wait and dynamic teams disabled.

| File | Records | Contents |
| --- | ---: | --- |
| ab.csv | 108 | Three paired trials per case, thread count and variant comparison |
| scaling.csv | 60 | Strong and weak scaling for V1 and V2 at five thread counts |
| counters.csv | 360 | Ten hardware events for each variant, shape, worker count and repetition |
| vectorization.csv | 3 | GCC diagnostics for the score-row loop |
| jobs.csv | 4 | Slurm state, hardware placement and elapsed time |
| topology.txt | | Mapping from sampled CPU IDs to sockets, NUMA domains and caches |
| environment.txt | | Hardware, compiler, perf version and OpenMP settings |

All successful jobs ran on the same AMD EPYC 9654 node. Baseline counters and the initial/V1 comparison used the same set of logical CPUs. Slurm assigned a different set to the V1/V2 comparison. Each A/B pair used the same cores in the same job. Job 838066 encountered a perf control-handshake error before recording a measurement. The fixed driver handled perf 6.9's NUL-terminated acknowledgement and job 838147 completed the initial profile.

The published CSVs concatenate the individual job files. Variant names and source digests identify the compiled translation units. The first two successful jobs labeled V1 as `candidate`; the published tables use `v1` throughout. The driver recorded maximum match gaps on the same deterministic inputs in the V2 run. Earlier records use the matching input digest to supply this field. The memory-model columns are calculated from the actual index representation and dimensions. The elapsed times, process RSS, input checksums, distances and counter values come from the runs. Each timed distance agreed with the two-row scalar oracle. Timings include one full warmed NvPD call. Counter values include all its threads in user space.

Each A/B speedup is the median of three time ratios paired by trial, input and thread count. On the square case, V1 ran 2.32× faster than the initial kernel on one core and 1.98× faster on eight. The 32-byte-vector loop retired 2.13 instructions per cell on one core, down from 21.76. V1 improved every measured shape and thread count, with median paired speedups from 1.45× to 5.85×. V2 reduced the steady four-letter index from 16 to 4 bytes per column but required a transient conversion buffer. Its speed ratios relative to V1 ranged from 0.88× to 1.14×.

In the fixed-size V1 scaling trials, eight cores took 27.9 ms and sixteen took 31.3 ms (three-trial medians). At sixteen cores the strong speedup was 5.26× over the one-worker NvPD implementation. With exactly 16,777,216 cells per worker, the weak-scaling time rose from 8.6 ms to 13.7 ms. The 4,194,304-column case is a separate cache diagnostic and exceeds the course dataset's maximum length.
