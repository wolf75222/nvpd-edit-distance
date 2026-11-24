#!/bin/bash -l
#SBATCH --job-name=chps-kokkos-cpu
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=kokkos-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false
srun --exact --mpi=none --cpu-bind=cores bash bench/kokkos/kokkos-step.sh cpu "results/kokkos/cpu-$SLURM_JOB_ID"
