#!/bin/bash -l
#SBATCH --job-name=nvpd-ctrace
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:55:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
bench=bench/cpu/chromosome-trace
out=results/trace/chromosomes/$SLURM_JOB_ID
if [[ -n ${SLURM_ARRAY_TASK_ID:-} ]]; then
    out="results/trace/chromosomes/$SLURM_ARRAY_JOB_ID/block-$SLURM_ARRAY_TASK_ID"
    export CHPS_TRACE_LONG_BLOCK="$SLURM_ARRAY_TASK_ID"
fi
mkdir -p "$out"
tmp=$(mktemp -d "$PWD/.trace.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
export TMPDIR="$tmp" OMP_DYNAMIC=false OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
cmake -S . -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native >"$tmp/config.log"
cmake --build "$tmp/build" --parallel 8 --target distances >"$tmp/build.log"
g++ -std=c++20 -O3 -march=native -fopenmp -Iinclude -Itests "$bench/trial.cpp" src/io/fasta.cpp "$tmp/build"/{libdistances,libedlib_static,libwfa_static}.a -lm -o "$tmp/trial"
{ hostname -s; g++ --version | head -1; lscpu | awk '/^Model name:|L1d cache:|L2 cache:|L3 cache:/'; } >"$out/environment.txt"
sha256sum src/cpu/split.cpp "$bench"/{run.sbatch,trial.cpp,measure.py} >"$out/sources.sha256"
(cd "$tmp" && sha256sum trial) >"$out/binary.sha256"
(cd "$1" && sha256sum chr1.fa chr18.fa chr20.fa) >"$out/inputs.sha256"
sha256sum "$2" | awk '{print $1,"chr21-2m.fa"}' >>"$out/inputs.sha256"
srun --cpu-bind=cores python3 "$bench/measure.py" "$tmp/trial" "$1" "$2" "$out/trials.csv"
