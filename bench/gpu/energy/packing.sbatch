#!/bin/bash -l
#SBATCH --job-name=nvpd-match-cache
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=12G
#SBATCH --time=00:17:00
#SBATCH --output=packing-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
set -u
inputs=$(realpath "${1:?provide the chromosome window directory}")
base=bench/gpu/energy
out=results/energy/gh200/$SLURM_JOB_ID
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-compact.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
flags=(-O3 -std=c++20 -arch=sm_90 -Iinclude -Isrc/cuda -Isrc/cpu -I"$base" -DCHPS_CUDA_TILE_ROWS=16)
nvcc "${flags[@]}" -c "$base/trial.cu" -o "$tmp/trial.o"
g++ -O3 -std=c++20 -fopenmp -Iinclude -c src/cpu/dp.cpp -o "$tmp/dp.o"
g++ -O3 -std=c++20 -Iinclude -c src/io/fasta.cpp -o "$tmp/fasta.o"
for kernel in direct packed; do
    packed=0
    [[ $kernel == direct ]] || packed=1
    nvcc "${flags[@]}" -DCHPS_CUDA_PACKED_MATCHES="$packed" src/cuda/cuda.cu "$tmp/"{trial,dp,fasta}.o -Xcompiler=-fopenmp,-pthread \
        -l:libnvidia-ml.so.1 -o "$tmp/$kernel"
done
{
    hostname -s
    nvidia-smi --query-gpu=name,uuid,pci.bus_id,driver_version,power.limit,memory.total --format=csv
    nvidia-smi --query-compute-apps=gpu_uuid,process_name --format=csv
    nvcc --version |tail -2
    sha256sum "$base"/{trial.cu,packing.sbatch,compare.py} src/cuda/{cuda.cu,local_tiles.cuh,packed_match.cuh}
} > "$out/environment.txt"
export OMP_NUM_THREADS=1
srun --mpi=none --cpu-bind=cores python3 "$base/compare.py" "$tmp" "$inputs" \
    results/kokkos/human/references.csv "$out"
