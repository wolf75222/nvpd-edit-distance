#!/bin/bash -l
#SBATCH --job-name=nvpd-schedules
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:50:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
spack load /lni3xp2
pmix_prefix=$(spack location -i /zp25yhu)
export LD_LIBRARY_PATH="$pmix_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set -u
out="results/mpi/schedules/$SLURM_JOB_ID"
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-schedules.XXXXXXXX")
shared=$(mktemp -d "$PWD/.schedules.XXXXXXXX")
trap 'rm -rf -- "$tmp" "$shared"' EXIT
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false CHPS_MPI_LAUNCHER=mpirun
CHPS_MPI_HOSTS=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export CHPS_MPI_HOSTS
cmake -S . -B "$tmp" -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native" >"$tmp/config.log"
cmake --build "$tmp" --target editdist-mpi --parallel 8 >"$tmp/build.log"
cp "$tmp/editdist-mpi" "$shared/editdist-mpi"
g++ -std=c++20 -O3 -DNDEBUG -march=native -fopenmp -Iinclude bench/mpi/schedules/oracles.cpp src/io/fasta.cpp \
    "$tmp/libdistances.a" "$tmp/libedlib_static.a" "$tmp/libwfa_static.a" -lm -o "$tmp/oracles"
sha256sum src/mpi/{mpi.cpp,scheduled_matrix.hpp} src/cpu/dp.cpp bench/mpi/schedules/{run.sbatch,measure.py,oracles.cpp,cohort.py} >"$out/sources.sha256"
(cd "$tmp" && sha256sum editdist-mpi oracles) >"$out/binaries.sha256"
{ printf 'nodes=%s\n' "$CHPS_MPI_HOSTS"; g++ --version | head -1; mpirun --version | head -1; } >"$out/environment.txt"
OMP_NUM_THREADS=32 "$tmp/oracles" data/mpi-schedules/locus.fa >"$out/scores.tsv"
python3 bench/mpi/schedules/measure.py "$shared/editdist-mpi" data/mpi-schedules "$out" "$shared" results/mpi/schedules/windows.tsv
