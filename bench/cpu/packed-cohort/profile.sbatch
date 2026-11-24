#!/bin/bash -l
#SBATCH --job-name=packed-pmu
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH --output=packed-pmu-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
root=$SLURM_SUBMIT_DIR
fasta=$(realpath "${1:?provide a cohort FASTA}")
out=$root/results/human/packed-pmu
perf=${PERF:-/project/r260105/rmdraux/tools/perf/bin/perf}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/packed-pmu.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$out"
cmake -S "$root" -B "$tmp/app" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native > "$tmp/configure.log"
cmake --build "$tmp/app" --target distances --parallel 2 > "$tmp/build.log"
g++ -std=c++20 -O3 -march=native -fopenmp -I"$root/include" "$root/bench/cpu/packed-cohort/profile.cpp" \
    "$root/src/cpu/packed.cpp" "$root/src/io/fasta.cpp" "$tmp/app"/{libdistances,libedlib_static,libwfa_static}.a -lm -o "$tmp/probe"
export LC_ALL=C OMP_NUM_THREADS=1
{
    echo "job=$SLURM_JOB_ID"
    hostname -s
    "$perf" --version
    lscpu | awk '/^Model name:/'
    (cd "$root" && sha256sum bench/cpu/packed-cohort/{profile.cpp,profile.sbatch} src/cpu/packed.cpp)
} > "$out/environment.txt"
mkfifo "$tmp/control" "$tmp/ack"
echo 'trial,method,pairs,hash,kernel_s,event,count,active_percent' > "$out/counters.csv"
methods=(byte packed edlib wfa2 one-sided)
for block in {0..5}; do
    for j in {0..4}; do
        method=${methods[$(( (block + j) % 5 ))]}
        srun -N1 -n1 -c1 --cpu-bind=cores "$perf" stat -x, -o "$tmp/counts" \
          --delay=-1 --control="fifo:$tmp/control,$tmp/ack" \
          -e '{cycles:u,instructions:u,L1-dcache-load-misses:u,branch-misses:u,rF045:u}' \
          "$tmp/probe" "$fasta" "$method" "$tmp/control" "$tmp/ack" > "$tmp/timing"
        grep -Eq '^[0-9]+,[0-9.e+-]+$' "$tmp/timing"
        awk -F, -v head="$block,$method,8128,$(cat "$tmp/timing")" \
          '$1 ~ /^[0-9]+$/ {print head "," $3 "," $1 "," $5}' "$tmp/counts" >> "$out/counters.csv"
    done
done
[ "$(wc -l < "$out/counters.csv")" -eq 151 ]
[ "$(cut -d, -f4 "$out/counters.csv" | tail -n+2 | sort -u | wc -l)" -eq 1 ]
