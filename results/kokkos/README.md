# Kokkos measurements on ROMEO

The Kokkos backend runs the NvPD cost recurrence as one range kernel per score row. The native CPU comparison uses the selected OpenMP SIMD kernel. The native GPU comparison uses the selected CUDA implementation with its tile optimization. The harness builds both compared backends in the same executable and tests the normal CLI separately. See the [build and benchmark procedure](../../bench/kokkos/kokkos.md).

`timings.csv` contains independent process trials for reproducible 5% synthetic inputs. Each trial warms the backend once, then records one complete cost-only call. All measured distances match edlib. `peak_rss_kib` is the process high-water mark. `device_bytes` counts the explicit score and index arrays on the accelerator. The CPU runs report zero device bytes.

On the EPYC 9654, Kokkos 4.5.01 selected OpenMP. The one-node CPU step of allocation 838222 used GCC 11.5 and eight assigned cores. Its Kokkos tests and timings completed before a subsequent independent MPI step failed. At 20,000 reference bases and eight workers, native OpenMP took 27.90 ms and Kokkos OpenMP took 110.22 ms (three-trial medians). The median paired native/Kokkos time ratio was 0.253. Seven additional CLI pairs checked cost and deterministic traceback, and a five-record matrix matched edlib for all ten unordered pairs.

On the GH200, Kokkos reported the Cuda execution space. Job 838230 used CUDA 12.6.3 and four host workers. At 20,000 reference bases, native CUDA took 17.06 ms and Kokkos CUDA took 44.27 ms. The median paired native/Kokkos time ratio was 0.385. At 5,000 bases, the corresponding medians were 4.16 ms and 10.75 ms. The native CTest and public CLI checks passed in both execution spaces.

The single cold CLI run in `human.json` compared 2,000,000 and 2,001,263 bases. Kokkos CUDA took 24.204 s and returned distance 5,599, matching the edlib and WFA2 references. It allocated 48,030,336 bytes of explicit device storage and reached 221,184 KiB of process peak RSS. The CLI created the Kokkos runtime before starting the timer.

The first accelerator build in 838220 encountered a CUDA 13 API incompatibility in Kokkos 4.5.01. Job 838227 built the library with CUDA 12.6.3, then passed edlib to nvcc_wrapper and hit an incompatible min/max definition. The final build used the Kokkos CMake launcher to route only dependent targets to nvcc_wrapper.

`builds.csv` records the launcher and step script digests used by each historical job. The public `bench/kokkos/kokkos.sbatch` is the portable equivalent of `recovery.sbatch`. It takes site-local paths through environment variables. `sources.sha256` identifies the CPU and GPU kernels used in the timings.

The native CUDA source in this comparison corresponds to commit `3f1840ff498bd68240cbf806def414ff32c1e331`, before the reusable matrix workspace. The isolated pair benchmark retains the same allocation scope.
