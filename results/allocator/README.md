# CUDA allocation reuse

The probe compares two allocation policies with the same production NvPD CUDA kernels. In normal mode every nonempty pair allocates and frees its match-index and score buffers. Reuse mode keeps one slot for each buffer over a serial batch. It retains a slot only when the next request has the same size, frees a mismatched slot before replacement, and releases both slots at batch exit. It never caches buffers for a list of sizes.

Run these commands from the project root on ROMEO. The checkout must contain commit `3f1840ff498bd68240cbf806def414ff32c1e331`. The natural cohort path is `data/phix174/cohort.fa`. The three fixed-size ACGT cohorts are generated in the job. Both jobs request one GH200 and four host workers. The toolchain is CUDA 13.0.2 with GPU target sm_90 and host OpenMP.

```bash
sbatch --account="$ROMEO_PROJECT" bench/allocator.sbatch
sbatch --account="$ROMEO_PROJECT" bench/allocator_validate.sbatch
```

The first launcher extracts the original source revision and checks the kernel digest. For a staged source tree, supply a Git archive with `ALLOCATOR_BASE_ARCHIVE`. Each job writes to `results/allocator/<job-id>/`. The second launcher builds the archived CLI and the current CLI, runs CTest on both, checks public traceback and tests a sequence of changes from wide to short rows and back to wide rows.

## Isolated probe

Job 838236 included the original CUDA translation unit in the probe. Only its calls to cudaMalloc and cudaFree were redirected through a counter. Both modes use the same instrumented binary and run the same serial alignments. Each process first warms one pair and releases its buffers. Batch timing includes the first allocation, all pair alignments and the final release. Independent edlib checks for every pair are outside the timer. Trials alternate the two modes, and each reported speedup is the median of three paired time ratios.

| Case | Pairs per batch | Physical cudaMalloc calls, normal → reuse | Peak reserved bytes, both modes | Median paired speedup |
| --- | ---: | ---: | ---: | ---: |
| Fixed 256 | 120 | 240 → 2 | 6,168 | 1.59× |
| Fixed 5,000 | 28 | 56 → 2 | 120,024 | 1.05× |
| Fixed 20,000 | 28 | 56 → 2 | 480,024 | 1.45× |
| ΦX174, 64 records | 2,016 | 4,032 → 250 | 150,892 | 1.38× |

`trials.csv` records the number of physical calls and reuse hits. Peak reserved bytes include idle retained slots. Peak live bytes count buffers held by a pair inside a call. The total allocated bytes sum all physical allocations in a batch. These measures cover the two explicit global-memory buffers and omit the CUDA context, driver-managed storage and per-kernel shared memory. Different shapes can leave idle storage between calls. The measurements establish equality of the two peaks only for the tested batches. Time attributed to cudaMalloc or cudaFree is host wall time inside these API calls and can include synchronization.

## Production CLI

A CUDA matrix owns one `chps::CudaWorkspace` and uses it for all serial pairs. A standalone `nvpd_cuda` call creates a private workspace and releases it on return. The row kernels and the choice of CUDA launch geometry are unchanged. The archived source and selected workspace source have the same SHA-256 digest for the block of code containing the kernels, `00bc56e4c8b5cf605776c53418175c33bebc3be023a5ca57f4606a5efc85275a`.
