#!/bin/bash -l
#SBATCH --job-name=cohort-ablation
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=ablation-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
fasta=$(realpath "${1:?provide the cohort FASTA}")
oracle=$(realpath "${2:?provide its edlib matrix}")
root=$(realpath "${3:-.}")
bench=$SLURM_SUBMIT_DIR/bench/cpu/surveillance
out=$SLURM_SUBMIT_DIR/results/phylogeny/sars-ablation
tmp=$(mktemp -d "${TMPDIR:-/tmp}/cohort-ablation.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$out"
cmake -S "$root" -B "$tmp/app" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native > "$tmp/configure.log"
cmake --build "$tmp/app" --target editdist --parallel 8 > "$tmp/build.log"
g++ -std=c++20 -O3 -march=native -fopenmp -I"$root/include" -I"$root/tests" "$bench/ablation.cpp" \
    "$root/src/io/fasta.cpp" "$tmp/app"/{libdistances,libedlib_static,libwfa_static}.a -lm -o "$tmp/ablation"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false
srun --mpi=none --ntasks=1 --cpu-bind=cores --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    "$tmp/app/editdist" matrix --input "$fasta" --algorithm nvpd-split --threads "$SLURM_CPUS_PER_TASK" > "$tmp/cli.tsv"
cmp "$tmp/cli.tsv" "$oracle"
srun --mpi=none --ntasks=1 --cpu-bind=cores --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    "$tmp/ablation" "$fasta" "$oracle" "$SLURM_CPUS_PER_TASK" > "$out/trials.csv"
{
    echo "job=$SLURM_JOB_ID"
    hostname -s
    lscpu | awk '/^Model name:|^Core\(s\) per socket:|^Socket\(s\):/'
    g++ --version | head -1
    sha256sum "$fasta" | awk '{print $1,"cohort.fa"}'
    sha256sum "$oracle" | awk '{print $1,"matrix.tsv"}'
    (cd "$tmp" && sha256sum ablation)
    (cd "$root" && sha256sum src/cpu/split.cpp include/chps/split.hpp src/cli/main.cpp src/libraries/wfa2/wfa_adapter.c)
    (cd "$bench" && sha256sum ablation.cpp run.sbatch)
} > "$out/environment.txt"
