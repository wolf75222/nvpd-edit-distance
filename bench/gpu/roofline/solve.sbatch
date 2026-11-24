#!/bin/bash -l
#SBATCH --job-name=nvpd-roof-solve
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=16G
#SBATCH --time=00:55:00
#SBATCH --output=roof-solve-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
spack load /phhwuh3
set -u
inputs=$(realpath "${1:?provide FASTA directory prepared by bench/kokkos/human/inputs.py}")
mode=${2:-fusion}
base=bench/gpu/roofline
out=${ROOF_OUT:-results/roofline/gpu/$SLURM_JOB_ID}
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/roof-solve.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
if [[ -f $base/baseline.cu ]]; then
    cp "$base/baseline.cu" "$tmp/baseline.cu"
else
    git show b5d08d1e5bb6d6b8e744bc6680b404caf569e991:src/cuda/cuda.cu > "$tmp/baseline.cu"
fi
python3 "$base/instrument.py" "$tmp/baseline.cu" "$tmp/experiment.cu"
sources=("$tmp/experiment.cu")
extra=()
measure=()
if [[ $mode == local || $mode == classes ]]; then
    cmake -S . -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90
    cmake --build "$tmp/build" --target editdist check_distances -j4
    OMP_NUM_THREADS=2 ctest --test-dir "$tmp/build" -R '^distances$' --output-on-failure |sed "s|$tmp|build|g" >"$out/ctest.txt"
    nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -Iinclude -Isrc/cpu -Isrc/cuda -I"$base" \
        -Dnvpd_cuda=experiment_cuda -DCudaWorkspace=ExperimentalWorkspace -c "$tmp/experiment.cu" -o "$tmp/experiment.o"
    sources=("$tmp/experiment.o")
    extra=(-DCHPS_ROOF_LOCAL)
    measure=(--local)
    if [[ $mode == classes ]]; then
        python3 "$base/instrument_classes.py" src/cuda "$tmp"
        nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -Iinclude -Isrc/cpu -I"$tmp" -I"$base" -Isrc/cuda \
            -Dnvpd_cuda=class_global -Dnvpd_cuda_local=class_local -DCudaWorkspace=ClassWorkspace \
            -c "$tmp/class.cu" -o "$tmp/class.o"
        sources+=("$tmp/class.o")
        extra+=(-DCHPS_ROOF_CLASSES)
        measure+=(--classes)
    fi
    python3 "$base/local_check.py" "$tmp/build/editdist" "$tmp" >"$out/cli.txt"
else
    cmake -S . -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
    cmake --build "$tmp/build" --target distances -j4
fi
nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -Iinclude -Isrc/cpu -Isrc/cuda -I"$base" \
    "${extra[@]}" "$base/trial.cu" "${sources[@]}" src/io/fasta.cpp "$tmp/build/"lib{distances,edlib_static,wfa_static}.a -lm -o "$tmp/trial"
"$tmp/trial" --check > "$out/validation.txt"
{
    hostname -s
    nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv
    nvcc --version |tail -2
    ncu --version |tail -1
    nsys --version
    grep -i ProfilingAdminOnly /proc/driver/nvidia/params || true
    sha256sum "$base"/{trial.cu,instrument.py,measure.py,solve.sbatch} "$base/paired_rows.cuh"
    if [[ $mode != fusion ]]; then sha256sum src/cuda/{cuda.cu,local_tiles.cuh} include/chps/cuda.hpp "$base/local_check.py"; fi
    if [[ $mode == classes ]]; then sha256sum "$base"/{width_classes.cuh,instrument_classes.py}; fi
    sha256sum "$tmp/experiment.cu" "$tmp/trial" |sed -E 's@  .*/@  @'
} > "$out/environment.txt"
export OMP_NUM_THREADS=1
srun --mpi=none --cpu-bind=cores python3 "$base/measure.py" "$tmp/trial" "$inputs" "$out" "$tmp" \
    results/kokkos/human/references.csv "${measure[@]}"
