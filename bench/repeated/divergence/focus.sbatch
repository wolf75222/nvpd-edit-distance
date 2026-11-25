#!/bin/bash -l
#SBATCH --job-name=nvpd-focus
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=00:58:00
#SBATCH --output=nvpd-focus-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
baseline=$(realpath "${1:?path to baseline repeat_kernel}")
out=${2:?new output directory}
mode=${3:?focus, low or select}
case $mode in
    focus) groups=(original-five held-out control-5pct);;
    low) groups=(100k-2pct 1m-2pct);;
    select) groups=(original-five held-out control-5pct 100k-2pct 1m-2pct);;
    *) exit 2;;
esac
selected_file=${CHPS_NVPD_SELECTION:-}
if [[ $mode == select && ${SLURM_RESTART_COUNT:-0} -gt 0 && -s "$out/sweep/selected-flags.txt" ]]; then
    selected_file=$(realpath "$out/sweep/selected-flags.txt")
    out="$out/resume-$SLURM_RESTART_COUNT"
fi
mkdir "$out"
out=$(realpath "$out")
work=$(mktemp -d /dev/shm/chps-focus.XXXXXXXX)
trap 'rm -rf -- "$work"' EXIT
build="$work/build"
bench="$PWD/bench/repeated/divergence"
initial=()
if [[ $mode != select ]]; then
    initial+=("-DCMAKE_CUDA_FLAGS=-DCHPS_NVPD_CORRIDOR_ROWS=32 -DCHPS_NVPD_CORRIDOR_COLUMNS=128 -DCHPS_NVPD_CORRIDOR_ANCHORS=0")
fi
cmake -S . -B "$build" "${initial[@]}" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DCHPS_CUDA=ON \
    -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_WARNINGS_AS_ERRORS=ON \
    -DCMAKE_PROJECT_INCLUDE="$PWD/bench/repeated/targets.cmake" >"$out/configure.log"
cmake --build "$build" --parallel 4 --target repeat_kernel editdist check_distances focus_cohort >"$out/build.log"
export OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
step=(srun --mpi=none --cpu-bind=cores)
if [[ $mode == select ]]; then
    python3 scripts/datasets/synthetic.py --out "$work/selection" --lengths 1000000 --rates 0.02 --pairs 11
    selection=("$work"/selection/pairs/*_p11_*.fa)
    control=data/synthetic/pairs/n100000_r0.1_p1_79f31facc7f5111c.fa
    score=$("$build/editdist" pair --algorithm edlib --threads 1 --input "${selection[0]}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["distance"])')
    if [[ -z $selected_file ]]; then
        "${step[@]}" bash "$bench/corridor_sweep.sh" "$build" "${selection[0]}" "$score" "$control" 9589 "$out/sweep"
        selected_file="$out/sweep/selected-flags.txt"
    fi
    cp "$selected_file" "$out/configuration.txt"
    read -r -a flags < "$selected_file"
    cmake -S . -B "$build" "-DCMAKE_CUDA_FLAGS:STRING=${flags[*]}" >> "$out/configure.log"
    cmake --build "$build" --parallel 4 --target repeat_kernel editdist check_distances focus_cohort >> "$out/build.log"
    c++ -std=c++20 -O2 -Iinclude "${flags[@]}" "$bench/corridor_plan.cpp" src/io/fasta.cpp -o "$work/plan"
    "$work/plan" "$control" "$work"/selection/pairs/*_p[1-5]_*.fa > "$out/geometry.csv"
fi
if [[ $mode == select ]] && command -v compute-sanitizer >/dev/null; then
    python3 scripts/datasets/synthetic.py --out "$work/memcheck" --lengths 32769 --rates 0.02 --pairs 1
    fasta=("$work"/memcheck/pairs/*.fa)
    score=$("$build/editdist" pair --algorithm edlib --threads 1 --input "${fasta[0]}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["distance"])')
    printf 'run\n' | "${step[@]}" compute-sanitizer --tool memcheck --error-exitcode 7 \
        "$build/repeat_kernel" "${fasta[0]}" nvpd-cuda "$score" > "$out/memcheck.txt" 2>&1
    # The short pair can fall below the cutoff after an edit changes its length.
    # Exercise the selected corridor separately on the 100 knt control.
    device_bytes=$(awk -F, 'NR == 2 && $4 == "corridor" { print $10 }' "$out/geometry.csv")
    test -n "$device_bytes"
    printf 'run\n' | "${step[@]}" compute-sanitizer --tool memcheck --error-exitcode 7 \
        "$build/repeat_kernel" "$control" nvpd-cuda 9589 > "$out/memcheck-corridor.txt" 2>&1
    awk -F, -v peak="$device_bytes" '$3 == 9589 && $6 == peak {matched=1} END {exit !matched}' \
        "$out/memcheck-corridor.txt"
fi
"$build/check_distances" > "$out/tests.txt"
{
    hostname -s
    printf 'allocation=%s\n' "$SLURM_JOB_ID"
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
    nvcc --version | tail -2
    g++ --version | head -1
    nsys --version
    sha256sum src/cuda/{cuda.cu,corridor.cuh,corridor_plan.hpp,corridor_anchor.hpp,local_tiles.cuh,packed_index.hpp,packed_match.cuh}
} > "$out/environment.txt"
"${step[@]}" python3 "$bench/focus.py" "$baseline" "$build/repeat_kernel" "$build/editdist" "$out/timing" \
    --groups "${groups[@]}" --blocks 8
if [[ $mode != focus ]]; then
    "${step[@]}" python3 "$bench/focus_memory.py" "$baseline" "$build/repeat_kernel" "$build/editdist" \
        "$out/memory" --groups 100k-2pct 1m-2pct --blocks 3
fi
if [[ $mode != low ]]; then
    "${step[@]}" python3 "$bench/cohort_memory.py" "$build/focus_cohort" "$build/editdist" "$out/cohort" --blocks 5
    "${step[@]}" python3 "$bench/activity.py" "$build/repeat_kernel" \
        data/synthetic/pairs/n100000_r0.1_p1_79f31facc7f5111c.fa 9589 "$out/activity.csv"
fi
