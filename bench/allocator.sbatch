#!/bin/bash -l
#SBATCH --job-name=chps-alloc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --constraint=armgpu
#SBATCH --time=00:10:00
#SBATCH --output=alloc-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
spack load /phhwuh3
set -u
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-alloc.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
out=results/allocator/$SLURM_JOB_ID
mkdir -p "$out"
base=$work/initial
mkdir -p "$base"
if [[ -n ${ALLOCATOR_BASE_ARCHIVE:-} ]]; then
    tar -xf "$ALLOCATOR_BASE_ARCHIVE" -C "$base"
else
    git archive 3f1840ff498bd68240cbf806def414ff32c1e331 CMakeLists.txt include src tests | tar -xC "$base"
fi
expected=059ab7bed0285b4d038c8edb40742877b0eed4386166eec8631238d36855950b
test "$(sha256sum "$base/src/cuda.cu" | cut -d' ' -f1)" = "$expected"
if [[ -d deps ]]; then ln -s "$PWD/deps" "$base/deps"; fi
cmake -S "$base" -B "$work/build" -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=OFF
cmake --build "$work/build" --target distances --parallel "$SLURM_CPUS_PER_TASK"
mkdir -p "$base/bench"
cp bench/allocator_probe.cu "$base/bench/"
g++ -std=c++20 -O3 -I"$base/include" -c "$base/src/fasta.cpp" -o "$work/fasta.o"
nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -I"$base/include" "$base/bench/allocator_probe.cu" \
    "$work/fasta.o" "$work/build/libdistances.a" "$work/build/libedlib_static.a" \
    "$work/build/libwfa_static.a" -lm -o "$work/probe"
(cd "$base" && sha256sum src/cuda.cu bench/allocator_probe.cu) > "$out/sources.sha256"
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv > "$out/device.txt"
nvcc --version >> "$out/device.txt"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
srun --mpi=none --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    python3 bench/allocator_measure.py "$work/probe" "$out/trials.csv"
