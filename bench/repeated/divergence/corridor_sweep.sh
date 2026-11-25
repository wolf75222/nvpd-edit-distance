#!/usr/bin/env bash
# Run inside one GH200 allocation with a disposable, configured project build.
set -euo pipefail
build=$(realpath "${1:?directory containing repeat_kernel}")
large=$(realpath "${2:?held-out 1M/2 percent FASTA}")
large_score=${3:?Edlib oracle}
control=$(realpath "${4:?historical 100k/10 percent FASTA}")
control_score=${5:?Edlib oracle}
mkdir "${6:?new output directory}"
out=$(realpath "$6")
root=$(pwd -P)
bench="$root/bench/repeated/divergence"
base_flags=$(sed -n 's/^CMAKE_CUDA_FLAGS:STRING=//p' "$build/CMakeCache.txt")
work=$(mktemp -d "${TMPDIR:-/dev/shm}/nvpd-sweep.XXXXXXXX")
restore_build() {
    local status=$?
    trap - EXIT
    if ! cmake -S "$root" -B "$build" "-DCMAKE_CUDA_FLAGS:STRING=$base_flags" >>"$out/build.log" ||
       ! cmake --build "$build" --target repeat_kernel --parallel 4 >>"$out/build.log"; then
        status=1
    fi
    rm -rf -- "$work"
    exit "$status"
}
trap restore_build EXIT
geometries=(32x128 16x64 16x128 32x64 32x256 64x128 64x256 16x256)
compile() {
    local geometry=$1 anchored=$2 name=$1
    if [[ $anchored == 1 ]]; then name=$geometry-anchor; fi
    local defs=(-DCHPS_NVPD_CORRIDOR_ROWS="${geometry%x*}" -DCHPS_NVPD_CORRIDOR_COLUMNS="${geometry#*x}" -DCHPS_NVPD_CORRIDOR_ANCHORS="$anchored")
    cmake -S "$root" -B "$build" "-DCMAKE_CUDA_FLAGS:STRING=$base_flags ${defs[*]}" >>"$out/build.log"
    cmake --build "$build" --target repeat_kernel --parallel 4 >>"$out/build.log"
    cp "$build/repeat_kernel" "$work/$name"
    c++ -O2 -std=c++20 -Iinclude "${defs[@]}" "$bench/corridor_plan.cpp" src/io/fasta.cpp -o "$work/plan"
    "$work/plan" "$large" "$control" >"$out/plan-$name.csv"
    (cd "$work" && sha256sum "$name") >>"$out/binaries.sha256"
}
for geometry in "${geometries[@]}"; do
    compile "$geometry" 0
done
export OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
python3 "$bench/corridor_sweep.py" "$work" "$out" "$large" "$large_score" "$control" "$control_score" "${geometries[@]}"
selected=$(cat "$out/selected.txt")
compile "$selected" 1
python3 "$bench/corridor_ablation.py" "$work" "$out" "$large" "$large_score" "$control" "$control_score" "$selected"
winner=$(cat "$out/selected-variant.txt")
for variant in $(printf '%s\n' 32x128 "$selected" "$winner" | sort -u); do
    python3 "$bench/activity.py" "$work/$variant" "$large" "$large_score" "$out/activity-$variant.csv" >"$out/nsys-$variant.log" 2>&1
done
