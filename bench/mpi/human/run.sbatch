#!/bin/bash -l
#SBATCH --job-name=chps-mpi-human
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:50:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
spack load /lni3xp2
pmix_prefix=$(spack location -i /zp25yhu)
export LD_LIBRARY_PATH="$pmix_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set -u
out="${CHPS_HUMAN_MPI_OUT:-results/mpi/human/$SLURM_JOB_ID}"
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/chps-human-mpi.XXXXXXXX")
shared=$(mktemp -d "$PWD/.humanmpi.XXXXXXXX")
trap 'rm -rf -- "$tmp" "$shared"' EXIT
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false CHPS_MPI_LAUNCHER=mpirun
CHPS_MPI_HOSTS=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export CHPS_MPI_HOSTS
cmake -S . -B "$tmp" -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native"
cmake --build "$tmp" --target editdist-mpi editdist --parallel 8
cp "$tmp/editdist-mpi" "$shared/editdist-mpi"
{ printf 'nodes=%s\n' "$CHPS_MPI_HOSTS"; g++ --version | head -1; mpirun --version | head -1; } >"$out/environment.txt"
sha256sum src/mpi/mpi.cpp src/cpu/dp.cpp bench/mpi/human/{run.sbatch,measure.py,windows.py} >"$out/sources.sha256"
python3 bench/mpi/human/measure.py "$shared/editdist-mpi" "$tmp/editdist" "${CHPS_HUMAN_MPI_INPUTS:-data/mpi-human}" "$out" "$shared"
