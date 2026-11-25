#!/bin/bash -l
#SBATCH --job-name=nvpd-launch
#SBATCH --partition=instant
#SBATCH --constraint=armgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --time=00:25:00
#SBATCH --output=nvpd-launch-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
base=$(realpath "${1:?baseline source at commit 9a9e23e or its repeat_kernel executable}")
mkdir "${2:?new output directory}"
out=$(realpath "$2")
work=$(mktemp -d /dev/shm/nvpd-launch.XXXXXXXX)
trap 'rm -rf -- "$work"' EXIT
opts=(-DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_WARNINGS_AS_ERRORS=ON)
if [[ -d $base ]]; then
    cmake -S "$base" -B "$work/base" "${opts[@]}" -DBUILD_TESTING=OFF \
        "-DCMAKE_PROJECT_INCLUDE=$base/bench/repeated/targets.cmake" >"$out/configure-base.log"
    cmake --build "$work/base" --parallel 4 --target repeat_kernel >"$out/build-base.log"
    baseline="$work/base/repeat_kernel"
else
    baseline=$base
    [[ -x $baseline ]]
fi
build="$work/candidate"
cmake -S . -B "$build" "${opts[@]}" -DBUILD_TESTING=ON \
    "-DCMAKE_PROJECT_INCLUDE=$PWD/bench/repeated/targets.cmake" >"$out/configure.log"
cmake --build "$build" --parallel 4 --target repeat_kernel editdist check_distances >"$out/build.log"
export OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
step=(srun --mpi=none --cpu-bind=cores)
bench="$PWD/bench/repeated/divergence"
{
    hostname -s
    printf 'job=%s\nbaseline_source_revision=9a9e23e\n' "$SLURM_JOB_ID"
    nvidia-smi --query-gpu=name,pci.bus_id,memory.total,driver_version,compute_mode --format=csv
    nvcc --version | tail -2
    g++ --version | head -1
    nsys --version
    sha256sum bench/repeated/kernel.cpp "bench/repeated/divergence/launch/measure.py"
} >"$out/environment.txt"
"${step[@]}" "$build/check_distances" >"$out/tests.txt"
"${step[@]}" python3 "$bench/launch/measure.py" "$baseline" "$build/repeat_kernel" "$build/editdist" "$out/timing"

python3 scripts/datasets/synthetic.py --out "$work/probes" --lengths 34000 100000 --rates 0.02 --pairs 1
short=("$work"/probes/pairs/n34000_*.fa)
long=("$work"/probes/pairs/n100000_*.fa)
oracle() {
    "$build/editdist" pair --input "$1" --algorithm edlib --threads 1 | python3 -c 'import json,sys; print(json.load(sys.stdin)["distance"])'
}
short_score=$(oracle "${short[0]}")
long_score=$(oracle "${long[0]}")
printf 'run\n' | "${step[@]}" compute-sanitizer --tool memcheck --error-exitcode=8 \
    "$build/repeat_kernel" "${short[0]}" nvpd-cuda-persistent "$short_score" >"$out/memcheck.txt" 2>&1
"${step[@]}" python3 "$bench/activity.py" "$baseline" "${long[0]}" "$long_score" "$out/activity-baseline.csv"
"${step[@]}" python3 "$bench/activity.py" "$build/repeat_kernel" "${long[0]}" "$long_score" \
    "$out/activity-persistent.csv" --algorithm nvpd-cuda-persistent
