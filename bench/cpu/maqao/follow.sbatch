#!/bin/bash -l
#SBATCH --job-name=chps-cqa-hints
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=x64cpu
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=maqao-%j.log
set -euo pipefail
base=$(realpath "${1:?pre-MAQAO source tree}")
prod=$(realpath "${2:?first production source tree}")
candidates=$(realpath "${3:?follow-up variant directory}")
inputs=$(realpath "${4:?prepared FASTA pairs}")
checks=$(realpath "${5:?directory containing check_variant.cpp}")
mkdir "${6:?new output directory}"
out=$(realpath "$6")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-cqa-hints.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=4
maqao=${MAQAO:-"$HOME/.local/bin/maqao"}
flags=(-O3 -g -fopenmp -std=c++20 -Wall -Wextra -Werror -Wpedantic -I"$prod/include")
names=(baseline aligned index32 unroll2 unroll4)
for mode in native generic; do
    mkdir "$work/$mode"
    arch=()
    [[ $mode == generic ]] || arch=(-march=native)
    for v in original "${names[@]}"; do
        src="$candidates/$v"
        [[ $v != original ]] || src="$base/src/cpu"
        [[ $v != baseline ]] || src="$prod/src/cpu"
        g++ "${arch[@]}" "${flags[@]}" -I"$src" "$src/dp.cpp" "$prod/src/io/fasta.cpp" \
            bench/cpu/maqao/kernel.cpp -o "$work/$mode/$v"
        [[ $v != original ]] || continue
        g++ "${arch[@]}" "${flags[@]}" -I"$src" "$src/dp.cpp" "$checks/check_variant.cpp" -o "$work/check"
        srun --mpi=none --cpu-bind=cores "$work/check" >> "$out/validation.txt"
    done
done
{
    hostname -s
    g++ --version | head -1
    "$maqao" --version | head -1
    echo "flags=${flags[*]}; native adds -march=native"
    srun --cpu-bind=cores bash -c 'grep Cpus_allowed_list /proc/self/status'
    sha256sum "$prod/src/cpu/nvpd_avx512.hpp" "$candidates"/*/nvpd_avx512.hpp
} > "$out/environment.txt"
time_calls() {
    srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$@" < /dev/null
}
profile() {
    local v=$1 threads=$2 iterations=$3
    local key="$v-$threads" xp="$work/xp-$v-$threads"
    [[ ! -e "$out/$key.tar.gz" ]] || return 0
    OMP_NUM_THREADS=$threads srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html \
        --lprof-params='sampling-rate=high' -xp="$xp" -- "$work/native/$v" "$inputs/near_32768.fa" \
        32768 "$threads" "$iterations" < /dev/null > "$out/oneview-$key.log" 2>&1
    tar -C "$xp" -czf "$out/$key.tar.gz" .
}
bins=()
for v in "${names[@]}"; do bins+=("$work/native/$v"); done
time_calls "$out/screen.csv" "${bins[@]}" --lengths 16384 32768 --blocks 10
python3 bench/cpu/maqao/select_variant.py "$out/screen.csv" > "$out/selection.csv"
mapfile -t finalists < <(tail -n +2 "$out/selection.csv" | cut -d, -f2 | sort -u)
bins=("$work/native/baseline" "$work/native/original")
generic=("$work/generic/baseline")
for v in "${finalists[@]}"; do
    [[ $v != baseline ]] || continue
    bins+=("$work/native/$v")
    generic+=("$work/generic/$v")
done
# With one finalist, each of the three rotated orders occurs three times.
time_calls "$out/long.csv" "${bins[@]}" --lengths 100000 --blocks 9
time_calls "$out/short.csv" "${bins[@]}" --lengths 32 128 512 1024 4096 --blocks 9 --calls 12
if (( ${#generic[@]} > 1 )); then
    time_calls "$out/generic.csv" "${generic[@]}" --lengths 32768 --blocks 9
fi
# Capture both four-core long profiles with the same number of calls.
profile original 4 437
profile baseline 4 437
# Each selected build is profiled on both thread counts; this also qualifies a shared implementation.
for v in "${finalists[@]}"; do
    profile "$v" 1 200
    profile "$v" 4 500
done
srun --mpi=none --cpu-bind=cores bash bench/cpu/maqao/lengths.sh \
    "$prod" "$candidates/aligned" "$inputs" "$out/lengths.csv"
