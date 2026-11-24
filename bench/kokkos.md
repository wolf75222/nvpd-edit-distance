# Kokkos NvPD

Build the application against an installed Kokkos 4 or later. The library configuration selects the default execution space. For a CUDA build of Kokkos, select its host C++ compiler when configuring the application. The Kokkos CMake launcher routes its dependent targets through `nvcc_wrapper`. This keeps ordinary host libraries on their C++ compiler.

```bash
cmake -S . -B build-kokkos -DCMAKE_BUILD_TYPE=Release \
    -DCHPS_KOKKOS=ON -DCMAKE_PREFIX_PATH="$KOKKOS_PREFIX"
cmake --build build-kokkos --parallel
build-kokkos/editdist pair --algorithm nvpd-kokkos --input pair.fa --threads 8
build-kokkos/editdist matrix --algorithm nvpd-kokkos --input cohort.fa --threads 8
```

The pair command reports the execution-space name. Thread count controls host workers. The backend builds the last-occurrence table on the host, copies it to a Kokkos View and submits one RangePolicy kernel per score row to the same execution-space instance. Final synchronization precedes the host result copy. The cost-only mode allocates two score rows. The traceback mode retains the score matrix for inputs of at most 20,000 symbols, copies it to the host and uses the shared CIGAR routine. The process RSS includes the host copy of the index. The reported device bytes count explicit score and index buffers in non-host Kokkos memory.

## Repeat the comparison

Start in the project root on an allocated Linux CPU or CUDA node. Set `KOKKOS_SOURCE` to a Git checkout of Kokkos that contains tag `4.5.01`. The script extracts upstream commit `175257a51ff29a0059ec48bcd233ee096b2c0438` into a temporary directory. Alternatively, set `KOKKOS_PREFIX` to an existing installation and set `CXX` to its compiler. The CPU comparison uses OpenMP. The GPU comparison builds for Hopper (compute capability 9.0) with CUDA and OpenMP. Use CUDA 12.6.3 with the pinned Kokkos source on ROMEO. Kokkos 4.5.01 still calls the CUDA 12 version of `cudaMemPrefetchAsync`; it failed to compile with CUDA 13. The script also sets the temporary wrapper’s probe architecture to Hopper.

```bash
bash bench/kokkos-step.sh cpu results/kokkos-cpu
bash bench/kokkos-step.sh gpu results/kokkos-gpu
```

Each run builds the normal CLI, runs the differential CTest suite and checks cost, traceback and matrix mode through the CLI. It generates the same 5,000 and 20,000-base pairs with the course generator at 5% input edit rate. Each backend warms up once in a fresh process, then times one complete distance call including allocation and data transfers. The two backends alternate order across three trials. CPU timings cover one thread and the allocated team size (up to eight). GPU trials use one device. The program checks all distances against edlib, records process peak RSS, and attributes the run to its reported execution space. Check `timings.csv`, `checks.json` and `environment.txt` in the output directory.

On ROMEO, the accelerator launcher loads the tested CUDA and CMake packages. To include the human case, pass the path to the pinned two-record FASTA and its independently verified distance.

```bash
export KOKKOS_SOURCE=/path/to/kokkos
export CHPS_HUMAN_PAIR=data/human/chr21-2m.fa CHPS_HUMAN_EXPECTED=5599
sbatch --account="$ROMEO_PROJECT" bench/kokkos.sbatch
```
