#!/bin/bash -l
#SBATCH --job-name=chps-mpi-repeat
#SBATCH --constraint=x64cpu
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
spack load /lni3xp2
pmix_prefix=$(spack location -i /zp25yhu)
export LD_LIBRARY_PATH="$pmix_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set -u
tmp=$(mktemp -d "${TMPDIR:-/tmp}/chps-pmpi.XXXXXXXX")
shared=$(mktemp -d "$PWD/.pmpi.XXXXXXXX")
trap 'rm -rf -- "$tmp" "$shared"' EXIT
base=bench/mpi/collectives
study=${1:-collectives}
case "$study" in collectives|scaling|pipeline) ;; *) exit 2 ;; esac
out=results/mpi/$study/$SLURM_JOB_ID
mkdir -p "$out"
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false CHPS_MPI_LAUNCHER=mpirun
CHPS_MPI_HOSTS=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export CHPS_MPI_HOSTS
cmake -S . -B "$tmp" -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON -DBUILD_TESTING=OFF \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native"
cmake --build "$tmp" --target editdist-mpi editdist --parallel 8
cp "$tmp/editdist-mpi" "$shared/"

{ hostname; echo "$CHPS_MPI_HOSTS"; mpirun --version | head -2; g++ --version | head -1; } > "$out/environment.txt"
if [[ $study == collectives ]]; then
    mpicc -shared -fPIC -O2 "$base/probe.c" -o "$shared/probe.so"
    python3 "$base/run.py" "$shared/editdist-mpi" "$tmp/editdist" "$shared/probe.so" "$out/trials.csv"
elif [[ $study == pipeline ]]; then
    mpicc -shared -fPIC -O2 "$base/probe.c" -o "$shared/probe.so"
    export CHPS_MPI_DATA_TMP="$shared"
    python3 bench/mpi/mpi_checks.py --binary "$shared/editdist-mpi" --oracle "$tmp/editdist" \
        --launcher "mpirun --host ${CHPS_MPI_HOSTS%%,*}:192 -n 4 --map-by ppr:4:node:PE=2 --bind-to core" > "$out/validation.txt"
    sha256sum src/mpi/mpi.cpp bench/mpi/pipeline/measure.py bench/mpi/mpi_checks.py bench/mpi/mpi_benchmark.py > "$out/source.sha256"
    python3 bench/mpi/pipeline/measure.py "$shared/editdist-mpi" "$tmp/editdist" "$shared/probe.so" "$out"
else
    python3 bench/mpi/scaling.py "$shared/editdist-mpi" "$tmp/editdist" data/phix174/cohort.fa "$out/trials.csv"
fi
