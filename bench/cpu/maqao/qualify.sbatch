#!/bin/bash -l
#SBATCH --job-name=chps-cqa-layout
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=x64cpu
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=maqao-%j.log
set -euo pipefail
original=$(realpath "${1:?pre-MAQAO source tree}")
first=$(realpath "${2:?first production source tree}")
layout=$(realpath "${3:?candidate CPU source directory}")
inputs=$(realpath "${4:?matrix qualification inputs}")
checks=$(realpath "${5:?differential check source directory}")
mkdir "${6:?new output directory}"
out=$(realpath "$6")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-cqa-layout.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=4
maqao=${MAQAO:-maqao}
flags=(-O3 -g -std=c++20 -fopenmp -Wall -Wextra -Wpedantic -Werror -I"$first/include")
for v in original baseline layout; do
    src="$first/src/cpu"
    [[ $v != original ]] || src="$original/src/cpu"
    [[ $v != layout ]] || src="$layout"
    g++ -march=native "${flags[@]}" "$src/dp.cpp" "$first/src/io/fasta.cpp" bench/cpu/maqao/kernel.cpp -o "$work/$v"
done
for mode in generic native; do
    arch=()
    [[ $mode != native ]] || arch=(-march=native)
    g++ "${arch[@]}" "${flags[@]}" "$layout/dp.cpp" "$checks/check_variant.cpp" -o "$work/check"
    srun --mpi=none --cpu-bind=cores "$work/check" > "$out/validation-$mode.txt"
done
{
    hostname -s
    awk -F: '/model name/{print $2;exit}' /proc/cpuinfo
    g++ --version | head -1
    "$maqao" --version | head -1
    echo "flags=${flags[*]}; -march=native"
    srun --cpu-bind=cores bash -c 'grep Cpus_allowed_list /proc/self/status'
    sha256sum "$original/src/cpu/dp.cpp" "$first/src/cpu/nvpd_avx512.hpp" "$layout/nvpd_avx512.hpp"
} > "$out/environment.txt"
time_calls() {
    srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$@" < /dev/null
}
# Both chromosome pairs appear in each shape; the rectangle-transpose arms reverse the same two FASTA records.
time_calls "$out/threshold.csv" "$work/baseline" "$work/layout" --lengths 39999 40000 40001 65519 65536 65553 95003 12345 90019 --blocks 8
time_calls "$out/primary.csv" "$work/original" "$work/baseline" "$work/layout" --lengths 50003 45077 --blocks 9
time_calls "$out/large-four.csv" "$work/original" "$work/baseline" "$work/layout" --lengths 131072 --threads 4 --blocks 9
for threads in 1 4; do
    repeats=40
    [[ $threads != 4 ]] || repeats=120
    for v in baseline layout; do
        xp="$work/xp-$v-$threads"
        OMP_NUM_THREADS=$threads srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html \
            --lprof-params='sampling-rate=high' -xp="$xp" -- "$work/$v" "$inputs/near_95003.fa" \
            95003 "$threads" "$repeats" < /dev/null > "$out/oneview-$v-$threads.log" 2>&1
        tar -C "$xp" -czf "$out/$v-$threads.tar.gz" .
    done
done
srun --mpi=none --cpu-bind=cores bash bench/cpu/maqao/layout_widths.sh "$first" "$layout" "$inputs" "$out/widths.csv" "${7:?index-pitch candidate CPU directory}" "$checks"
