#!/bin/bash -l
#SBATCH --job-name=chps-edit-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
exec bash bench/romeo-run.sh cpu "$@"
