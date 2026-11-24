#!/bin/bash -l
#SBATCH --job-name=nvpd-candidate-loci
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=16G
#SBATCH --time=00:22:00
#SBATCH --output=loci-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
spack load /phhwuh3
set -u
inputs=$(realpath "${1:?provide chromosome window directory}")
out=results/energy/gh200/$SLURM_JOB_ID
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-loci.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
cmake -S . -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON \
    -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_CUDA_TILE_ROWS=16 -DCHPS_CUDA_PACKED_MATCHES=ON
cmake --build "$tmp/build" --target editdist check_distances -j4
export OMP_NUM_THREADS=1
ctest --test-dir "$tmp/build" -R '^distances$' --output-on-failure |sed "s|$tmp|build|" >"$out/ctest.txt"
python3 bench/gpu/roofline/local_check.py "$tmp/build/editdist" "$tmp" > "$out/cli-check.txt"
{
    hostname -s
    nvidia-smi --query-gpu=name,uuid,pci.bus_id,driver_version --format=csv
    nvcc --version |tail -2
    sha256sum src/cuda/{cuda.cu,local_tiles.cuh,packed_match.cuh} CMakeLists.txt bench/gpu/energy/loci.{py,sbatch}
} > "$out/environment.txt"
srun --mpi=none --cpu-bind=cores python3 bench/gpu/energy/loci.py "$tmp/build/editdist" \
    "$inputs" results/kokkos/human/references.csv "$out" "$tmp"
