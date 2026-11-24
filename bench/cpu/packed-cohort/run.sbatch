#!/bin/bash -l
#SBATCH --job-name=packed-cohort
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --mem=12G
#SBATCH --time=00:59:00
#SBATCH --output=packed-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
root=$SLURM_SUBMIT_DIR
fasta=$(realpath "${1:?provide a cohort FASTA}")
out=${CHPS_COHORT_RESULTS:-$root/results/human/packed-cohort}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/packed-cohort.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$out"
cmake -S "$root" -B "$tmp/app" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native > "$tmp/configure.log"
cmake --build "$tmp/app" --target distances --parallel 8 > "$tmp/build.log"
g++ -std=c++20 -O3 -march=native -fopenmp -I"$root/include" \
    "$root/bench/cpu/packed-cohort/trial.cpp" "$root/src/cpu/packed.cpp" "$root/src/io/fasta.cpp" \
    "$tmp/app"/{libdistances,libedlib_static,libwfa_static}.a -lm -o "$tmp/trial"
"$tmp/trial" --check
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false CHPS_ORACLE_CACHE="$tmp/oracle.bin"
{
    echo "job=$SLURM_JOB_ID"
    hostname -s
    lscpu | awk '/^Model name:|^Core\(s\) per socket:|^Socket\(s\):|^Thread\(s\) per core:/'
    g++ --version | head -1
    sha256sum "$fasta" | awk '{print $1,"cohort.fa"}'
    (cd "$root" && sha256sum bench/cpu/packed-cohort/{prepare.py,trial.cpp,run.sbatch} src/cpu/packed.cpp include/chps/packed.hpp)
} > "$out/environment.txt"
for nt in 64 16 1 4 32 8 2; do
    mode=full
    if [[ $nt != 64 ]]; then mode=scale; fi
    srun --mpi=none --ntasks=1 --cpu-bind=cores --cpus-per-task="$nt" "$tmp/trial" "$fasta" \
        "$out" "$nt" 5 "$mode" > "$out/trials-$nt.csv"
    if [[ $nt != 64 ]]; then
        srun --mpi=none --ntasks=1 --cpu-bind=cores --cpus-per-task="$nt" "$tmp/trial" "$fasta" \
            "$out" "$nt" 5 weak > "$out/weak-$nt.csv"
    fi
done
