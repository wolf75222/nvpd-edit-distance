#!/bin/bash -l
#SBATCH --job-name=nvpd-mpi-kokkos
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
spack load /lni3xp2
pmix=$(spack location -i /zp25yhu)
export LD_LIBRARY_PATH="$pmix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set -u
pin=175257a51ff29a0059ec48bcd233ee096b2c0438
tmp=$(mktemp -d "$PWD/.mpiworkers.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
export TMPDIR="$tmp" OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
export CHPS_MPI_LAUNCHER=mpirun
CHPS_MPI_HOSTS=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export CHPS_MPI_HOSTS
out="$PWD/results/mpi/workers/$SLURM_JOB_ID"
mkdir -p "$tmp/src" "$out"
git -C "${KOKKOS_SOURCE:?set KOKKOS_SOURCE to the upstream checkout}" archive "$pin" | tar -xf - -C "$tmp/src"
cmake -S "$tmp/src" -B "$tmp/kbuild" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 -DKokkos_ENABLE_OPENMP=ON -DCMAKE_INSTALL_PREFIX="$tmp/install" >"$tmp/kconfig.log"
cmake --build "$tmp/kbuild" -j 16 --target install >"$tmp/kbuild.log"
cmake -S . -B "$tmp/app" -DCMAKE_PREFIX_PATH="$tmp/install" -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON -DCHPS_KOKKOS=ON -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS=-march=native >"$tmp/config.log"
cmake --build "$tmp/app" -j 16 --target editdist-mpi >"$tmp/build.log"
sha256sum src/mpi/{mpi.cpp,scheduled_matrix.hpp} src/cpu/dp.cpp bench/mpi/workers/{compare.py,run.sbatch} >"$out/sources.sha256"
(cd "$tmp/app" && sha256sum editdist-mpi) >"$out/binary.sha256"
{ echo "nodes=$CHPS_MPI_HOSTS"; echo "kokkos=$pin"; g++ --version | head -1; mpirun --version | head -1; } >"$out/environment.txt"
(cd "$(dirname "$1")" && sha256sum "$(basename "$1")") >"$out/input.sha256"
(cd "$(dirname "$2")" && sha256sum "$(basename "$2")") >"$out/reference.sha256"
python3 bench/mpi/workers/compare.py "$tmp/app/editdist-mpi" "$1" "$2" "$out" "$tmp"
