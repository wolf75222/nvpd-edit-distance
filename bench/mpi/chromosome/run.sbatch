#!/bin/bash -l
#SBATCH --job-name=nvpd-chr1-matrix
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
spack load /lni3xp2
pmix=$(spack location -i /zp25yhu)
export LD_LIBRARY_PATH="$pmix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set -u
work=$(mktemp -d "$PWD/.chr1-matrix.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
out="$PWD/results/mpi/chromosome/$SLURM_JOB_ID"
mkdir -p "$out"
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false OMP_NUM_THREADS=8
export CHPS_MPI_LAUNCHER=mpirun
CHPS_MPI_HOSTS=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export CHPS_MPI_HOSTS
cmake -S . -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS=-march=native >"$work/config.log"
cmake --build "$work/build" --target editdist-mpi -j 16 >"$work/build.log"
sha256sum src/mpi/{mpi.cpp,scheduled_matrix.hpp,indexed_rows.hpp} src/cpu/dp.cpp bench/mpi/chromosome/{cohort.py,measure.py,run.sbatch} >"$out/sources.sha256"
(cd "$1" && sha256sum chr1.fa manifest.tsv) >"$out/inputs.sha256"
(cd "$work/build" && sha256sum editdist-mpi) >"$out/binary.sha256"
{ echo "nodes=$CHPS_MPI_HOSTS"; g++ --version | sed -n '1p'; mpirun --version | sed -n '1p'; lscpu | sed -n '1,30p'; } >"$out/environment.txt"
python3 bench/mpi/chromosome/measure.py "$work/build/editdist-mpi" "$1" "$out" "$work"
