#!/bin/bash -l
#SBATCH --job-name=chps-alloc-cli
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --constraint=armgpu
#SBATCH --time=00:10:00
#SBATCH --output=alloc-cli-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
spack load /phhwuh3
set -u
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-alloc-cli.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
out=results/allocator/$SLURM_JOB_ID
mkdir -p "$out"
baseline=$work/initial
mkdir "$baseline"
if [[ -n ${ALLOCATOR_BASE_ARCHIVE:-} ]]; then
    tar -xf "$ALLOCATOR_BASE_ARCHIVE" -C "$baseline"
else
    git archive 3f1840ff498bd68240cbf806def414ff32c1e331 CMakeLists.txt include src tests | tar -xC "$baseline"
fi
expected=059ab7bed0285b4d038c8edb40742877b0eed4386166eec8631238d36855950b
test "$(sha256sum "$baseline/src/cuda.cu" | cut -d' ' -f1)" = "$expected"
if [[ -d deps ]]; then ln -s "$PWD/deps" "$baseline/deps"; fi
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
for label in initial current; do
    root=$PWD
    if [[ $label == initial ]]; then root=$baseline; fi
    cmake -S "$root" -B "$work/$label-build" -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON \
        -DCMAKE_CUDA_ARCHITECTURES=90 -DBUILD_TESTING=ON
    cmake --build "$work/$label-build" --parallel "$SLURM_CPUS_PER_TASK"
    ctest --test-dir "$work/$label-build" --output-on-failure
    if [[ $label == initial ]]; then sources=(src/cuda.cu src/main.cpp); else sources=(src/cuda/cuda.cu src/cli/main.cpp); fi
    (cd "$root" && sha256sum "${sources[@]}" CMakeLists.txt) > "$out/$label.sha256"
done
srun --mpi=none --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    python3 bench/gpu/allocation/allocator_cli.py "$work/initial-build/editdist" "$work/current-build/editdist" "$out/production.csv"
