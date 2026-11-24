#!/bin/bash -l
#SBATCH --job-name=chps-hg002
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=12G
#SBATCH --time=00:25:00
#SBATCH --output=hg002-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
data=$(realpath "$1")
out=$(realpath -m "${2:-results/human/chromosome-dispatch/$SLURM_JOB_ID}")
mkdir -p "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-hg002.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
study=bench/cpu/chromosome-dispatch
cmake -S . -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native \
    -DCMAKE_PROJECT_INCLUDE="$PWD/$study/targets.cmake"
cmake --build "$work/build" --parallel 8 --target chromosome_trial
{ hostname -s; g++ --version | head -1; lscpu | awk '/Model name:|L1d cache:|L2 cache:|L3 cache:|NUMA node\(s\):/'; sha256sum src/cpu/{sparse_frontier,split}.cpp "$study"/{measure.py,trial.cpp}; } > "$out/environment.txt"
cp "$data/manifest.tsv" "$out/manifest.tsv"
srun --cpu-bind=cores python3 "$study/measure.py" "$work/build/chromosome_trial" "$data" "$out/trials.csv"
