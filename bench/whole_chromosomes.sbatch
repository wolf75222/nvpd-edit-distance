#!/bin/bash -l
#SBATCH --job-name=chps-whole-human
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
out=results/human/whole/$SLURM_JOB_ID
mkdir -p "$out"
cmake -S . -B build-whole -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build-whole --target editdist --parallel 1
export OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
hostname > "$out/environment.txt"
lscpu | awk -F: '/^Model name:|^Architecture:/{print}' >> "$out/environment.txt"
g++ --version | head -1 >> "$out/environment.txt"
sha256sum build-whole/editdist >> "$out/environment.txt"
srun --mpi=none --cpu-bind=cores python3 bench/whole_chromosomes.py build-whole/editdist results/human/whole/manifest.tsv "$out/trials.csv"
