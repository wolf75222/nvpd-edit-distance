#!/bin/bash -l
#SBATCH --job-name=chps-kokkos-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --constraint=armgpu
#SBATCH --time=00:15:00
#SBATCH --output=kokkos-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /uxlsycg
set -u
nvcc --version
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
srun --mpi=none --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    bash bench/kokkos/kokkos-step.sh gpu "results/kokkos/gpu-$SLURM_JOB_ID"
