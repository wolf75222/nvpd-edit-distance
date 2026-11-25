#!/bin/bash
# Run inside one four-core Slurm step to compare row strides between the two screen lengths.
set -euo pipefail
root=$(realpath "${1:?first production source tree}")
aligned=$(realpath "${2:?aligned candidate CPU directory}")
inputs=$(realpath "${3:?prepared input FASTAs}")
out=$(realpath -m "${4:?new CSV path}")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-cqa-lengths.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=4
flags=(-O3 -g -march=native -std=c++20 -fopenmp -Wall -Wextra -Wpedantic -Werror -I"$root/include")
for v in baseline aligned; do
    src="$root/src/cpu"
    [[ $v == baseline ]] || src="$aligned"
    g++ "${flags[@]}" "$src/dp.cpp" "$root/src/io/fasta.cpp" bench/cpu/maqao/kernel.cpp -o "$work/$v"
done
python3 bench/cpu/maqao/compare.py "$inputs" "$out" "$work/baseline" "$work/aligned" \
    --lengths 24001 40003 60007 80009 --blocks 8
