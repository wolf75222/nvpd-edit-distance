#!/bin/bash -l
#SBATCH --job-name=nvpd-energy-window
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --time=00:07:00
#SBATCH --output=sensor-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
set -u
inputs=$(realpath "${1:?provide chromosome window directory}")
out=results/energy/gh200/$SLURM_JOB_ID
base=bench/gpu/energy
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-sensor.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
flags=(-O3 -std=c++20 -arch=sm_90 -DCHPS_CUDA_PACKED_MATCHES=0 -Iinclude -Isrc/cuda -Isrc/cpu -DCHPS_CUDA_TILE_ROWS=16 -Xcompiler=-fopenmp,-pthread)
nvcc "${flags[@]}" "$base/geometry.cu" src/cpu/dp.cpp src/io/fasta.cpp -o "$tmp/geometry"
nvcc "${flags[@]}" "$base/trial.cu" src/cuda/cuda.cu src/cpu/dp.cpp src/io/fasta.cpp -l:libnvidia-ml.so.1 -o "$tmp/trial"
{
    hostname -s
    nvidia-smi --query-gpu=name,uuid,pci.bus_id,driver_version,power.limit --format=csv
    sha256sum "$base"/{trial.cu,geometry.cu,sensor.sbatch,sensor.py} src/cuda/{cuda.cu,local_tiles.cuh,packed_match.cuh}
} > "$out/environment.txt"
export OMP_NUM_THREADS=1
srun --mpi=none --cpu-bind=cores python3 "$base/sensor.py" "$tmp" "$inputs" "$out"
