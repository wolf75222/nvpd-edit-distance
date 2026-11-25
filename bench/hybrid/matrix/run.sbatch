#!/bin/bash -l
#SBATCH --job-name=chps-matrix
#SBATCH --constraint=armgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=matrix-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
inputs=${1:?directory containing FASTA inputs}
panel=${2:-human}
out="results/hybrid/jobs/$SLURM_JOB_ID"
mkdir -p "$out"
temp=$(mktemp -d "${TMPDIR:-/tmp}/chps-matrix.XXXXXXXX")
trap 'rm -rf -- "$temp"' EXIT
cmake -S . -B "$temp/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_CUDA_TILE_ROWS=16 \
    -DCHPS_CUDA_PACKED_MATCHES=ON -DCHPS_WARNINGS_AS_ERRORS=ON
cmake --build "$temp/build" --parallel "$SLURM_CPUS_PER_TASK" --target editdist-route
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=FALSE
if [[ $panel == human ]]; then cp "$inputs/manifest.tsv" "$out/manifest.tsv"; fi
{
    hostname -s
    g++ --version | head -1
    cmake --version | head -1
    nvcc --version | tail -2
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
    lscpu | awk '/^(Architecture|Model name):/'
    if [[ $panel == human ]]; then sha256sum "$inputs"/{near,mixed,distant}.fa
    else sha256sum "$inputs/surveillance.fa"; fi
    sha256sum bench/hybrid/matrix/{measure.py,run.sbatch} src/hybrid/* src/cuda/cuda.cu
    env | sort | sed -n '/^OMP_/p'
} > "$out/environment.txt"
srun --mpi=none --cpu-bind=cores python3 bench/hybrid/matrix/measure.py \
    "$temp/build/editdist-route" "$inputs" "$out" "$temp/matrices" --threads "$SLURM_CPUS_PER_TASK" --panel "$panel"
