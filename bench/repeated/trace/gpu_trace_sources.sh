#!/bin/bash
# Reconstruct the exact historical 5k experiment, including host and harness.
set -euo pipefail
out=${1:?expected destination}
mkdir -p "$out"
harness=ba186247f55f8f135499270450fddd0582d54ba3
git archive 4724f68bf7dbfb83b133422541b165b5c7124da1 src/dp.cpp src/dp_internal.hpp \
    src/fasta.cpp include/chps/dp.hpp include/chps/fasta.hpp | tar -x -C "$out"
git archive "$harness" bench/profile_gpu.cu | tar -x -C "$out"
git show f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83:src/cuda.cu > "$out/baseline.cu"
git show 21d7b7751f266ee60ce4ce44b1419d69963571db:src/cuda.cu > "$out/candidate.cu"
git show "$harness:results/profile/gpu/harness.sha256" > "$out/sources.sha256"
git show "$harness:results/profile/gpu/838170/variants.sha256" | sed 's/$/.cu/' >> "$out/sources.sha256"
(cd "$out" && sha256sum -c sources.sha256)
