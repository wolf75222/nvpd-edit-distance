# Kokkos measurements on ROMEO

The Kokkos backend runs the NvPD cost recurrence as one range kernel per score row. The native CPU comparison uses the selected OpenMP SIMD kernel. The native GPU comparison uses the selected CUDA implementation with its tile optimization. The harness builds both compared backends in the same executable and tests the normal CLI separately. See the [build and benchmark procedure](../../bench/kokkos.md).

`timings.csv` contains independent process trials for reproducible 5% synthetic inputs. Each trial warms the backend once, then records one complete cost-only call. All measured distances match edlib. `peak_rss_kib` is the process high-water mark. `device_bytes` counts the explicit score and index arrays on the accelerator. The CPU runs report zero device bytes.

On the EPYC 9654, Kokkos 4.5.01 selected OpenMP. The one-node CPU step of allocation 838222 used GCC 11.5 and eight assigned cores. Its Kokkos tests and timings completed before a subsequent independent MPI step failed. At 20,000 reference bases and eight workers, native OpenMP took 27.90 ms and Kokkos OpenMP took 110.22 ms (three-trial medians). The median paired native/Kokkos time ratio was 0.253. Seven additional CLI pairs checked cost and deterministic traceback, and a five-record matrix matched edlib for all ten unordered pairs.
