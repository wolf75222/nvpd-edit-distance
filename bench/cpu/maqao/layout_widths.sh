#!/bin/bash
# Keep each binary on the same four cores; compare neighboring logical widths on one EPYC.
set -euo pipefail
first=$(realpath "${1:?first production source tree}")
layout=$(realpath "${2:?candidate CPU directory}")
inputs=$(realpath "${3:?qualification pairs}")
out=$(realpath -m "${4:?new CSV}")
index_pitch=$(realpath "${5:-$(dirname "$layout")/index-pitch}")
check_dir=$(realpath "${6:-$(dirname "$layout")/variants}")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-cqa-widths.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=4
for v in baseline layout index-pitch; do
    src="$first/src/cpu"
    [[ $v != layout ]] || src="$layout"
    [[ $v != index-pitch ]] || src="$index_pitch"
    g++ -O3 -g -march=native -std=c++20 -fopenmp -Wall -Wextra -Wpedantic -Werror -I"$first/include" \
        "$src/dp.cpp" "$first/src/io/fasta.cpp" bench/cpu/maqao/kernel.cpp -o "$work/$v"
done
for mode in generic native; do
    arch=()
    [[ $mode == generic ]] || arch=(-march=native)
    g++ -O3 -std=c++20 -fopenmp -Wall -Wextra -Werror "${arch[@]}" -I"$first/include" \
        "$index_pitch/dp.cpp" "$check_dir/check_variant.cpp" -o "$work/check"
    "$work/check" > "$(dirname "$out")/validation-index-pitch-$mode.txt"
done
python3 bench/cpu/maqao/compare.py "$inputs" "$out" "$work/baseline" "$work/layout" "$work/index-pitch" \
    --lengths 49152 65521 65535 65536 65537 95003 98304 --blocks 9
