#!/bin/bash -l
#SBATCH --job-name=nvpd-refine
#SBATCH --partition=instant
#SBATCH --constraint=armgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --time=00:25:00
#SBATCH --output=nvpd-refine-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
reference=$(realpath "${1:?source tree at 148ecdd}")
candidate=$PWD
mkdir -p "${2:?output directory}"
out=$(realpath "$2")
work=$(mktemp -d /dev/shm/nvpd-refine.XXXXXXXX)
trap 'rm -rf -- "$work"' EXIT
export TMPDIR="$work" OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close OMP_WAIT_POLICY=ACTIVE
opts=(-DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_WARNINGS_AS_ERRORS=ON -DBUILD_TESTING=ON)
for archives in "${CHPS_DEPS:-$candidate/deps}" "$reference/deps"; do
    if [[ -f "$archives/edlib.tar.gz" && -f "$archives/wfa2.tar.gz" ]]; then break; fi
done
[[ -f "$archives/edlib.tar.gz" && -f "$archives/wfa2.tar.gz" ]] || { echo 'Missing pinned release archives; supply CHPS_DEPS' >&2; exit 2; }
archives=$(realpath "$archives")
for name in reference candidate; do
    tree=${!name}
    if [[ ! -e "$tree/deps" && ! -L "$tree/deps" ]]; then ln -s "$archives" "$tree/deps"; fi
    cmake -S "$tree" -B "$work/$name" "${opts[@]}" \
        "-DCMAKE_PROJECT_INCLUDE=$tree/bench/repeated/targets.cmake" >"$out/configure-$name.log"
    cmake --build "$work/$name" --parallel 4 --target repeat_kernel editdist check_hybrid_pair >"$out/build-$name.log" 2>&1 || {
        cat "$out/build-$name.log" >&2
        exit 2
    }
done
{
    hostname -s
    echo "job=$SLURM_JOB_ID"
    echo 'reference_source_revision=148ecdd'
    nvcc --version | tail -2
    g++ --version | head -1
    nvidia-smi --query-gpu=name,pci.bus_id,driver_version --format=csv
    sha256sum "$archives/edlib.tar.gz" "$archives/wfa2.tar.gz"
    for name in reference candidate; do
        (cd "${!name}" && sha256sum src/hybrid/pair.cpp src/cuda/cuda.cu bench/repeated/kernel.cpp) | sed "s|  |  $name/|"
    done
    sha256sum bench/repeated/divergence/launch/{measure,refine}.py
} >"$out/environment.txt"
step=(srun --mpi=none --cpu-bind=cores)
if [[ -n ${SLURM_STEP_ID:-} && ${SLURM_STEP_ID:-} != batch && ${SLURM_STEP_ID:-} != extern ]]; then step=(); fi
for name in reference candidate; do
    "${step[@]}" ctest --test-dir "$work/$name" -R '^hybrid_pair(_gpu)?$' --no-tests=error --output-on-failure >"$out/tests-$name.txt"
done
"${step[@]}" python3 bench/repeated/divergence/launch/refine.py run \
    "$work/reference" "$work/candidate" "$out/timing"
python3 bench/repeated/divergence/launch/refine.py analyze "$out/timing" "$out/analysis"
