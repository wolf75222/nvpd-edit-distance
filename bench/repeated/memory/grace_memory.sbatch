#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=grace-rss
#SBATCH --output=memory-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
build=$(mktemp -d "${TMPDIR:-/tmp}/grace-rss.XXXXXX")
trap 'rm -rf -- "$build"' EXIT
out=results/repeated/jobs/$SLURM_JOB_ID
mkdir -p "$out"
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 \
    "-DCMAKE_PROJECT_INCLUDE=$PWD/bench/repeated/targets.cmake"
cmake --build "$build" --target repeat_kernel --parallel 4
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
{
  hostname -s
  uname -sr
  getconf PAGESIZE
  g++ --version | head -1
  sha256sum bench/repeated/{inputs.csv,kernel.cpp} bench/repeated/memory/{grace_memory.py,grace_memory.sbatch}
} > "$out/environment.txt"
srun --mpi=none --cpu-bind=cores python3 bench/repeated/memory/grace_memory.py "$build/repeat_kernel" \
    results/repeated/gh200/summary.csv "$out"
