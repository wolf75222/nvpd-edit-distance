#!/bin/bash -l
#SBATCH --job-name=wfa-trace-ab
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=trace-ab-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
fasta=$(realpath "${1:?provide the cohort FASTA}")
oracle=$(realpath "${2:?provide the edlib matrix}")
out=$PWD/results/alternatives/wfa-portability/$SLURM_JOB_ID
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/wfa-trace-ab.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
cmake -S . -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native > "$tmp/config.log"
cmake --build "$tmp/build" --target distances --parallel 8 > "$tmp/build.log"
lib=$tmp/build
upstream=$lib/_deps/wfa2-src
cp "$lib/libwfa_static.a" "$tmp/previous.a"
units=(wavefront_backtrace wavefront_pcigar wavefront_extend_kernels_avx)
for unit in "${units[@]}"; do
    gcc -O3 -DNDEBUG -march=native -I"$upstream" -I"$upstream/wavefront" -I"$upstream/system" \
        -c "$upstream/wavefront/$unit.c" -o "$tmp/$unit.c.o"
    ar rs "$tmp/previous.a" "$tmp/$unit.c.o"
done
g++ -std=c++20 -O3 -DNDEBUG -march=native -fopenmp -Iinclude -Itests -c \
    bench/cpu/alternatives/wfa_trace_batch.cpp -o "$tmp/driver.o"
g++ -std=c++20 -O3 -DNDEBUG -march=native -fopenmp -Iinclude -c src/io/fasta.cpp -o "$tmp/fasta.o"
for variant in previous final; do
    archive=$lib/libwfa_static.a
    if [[ $variant == previous ]]; then archive=$tmp/previous.a; fi
    g++ -fopenmp "$tmp/driver.o" "$tmp/fasta.o" "$lib/libdistances.a" "$lib/libedlib_static.a" \
        "$archive" -lm -o "$tmp/$variant"
done
export OMP_DYNAMIC=false OMP_PLACES=cores OMP_PROC_BIND=close OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun --mpi=none --cpu-bind=cores bash -eu -c '
    tmp=$1 fasta=$2 oracle=$3
    echo block,position,mode,variant,seconds
    for block in {0..8}; do
        for mode in score trace; do
            variants=(previous final)
            if (( block % 2 )); then variants=(final previous); fi
            pos=0
            for variant in "${variants[@]}"; do
                pos=$((pos+1))
                "$tmp/$variant" "$fasta" "$mode" > "$tmp/matrix.tsv" 2> "$tmp/seconds"
                cmp "$oracle" "$tmp/matrix.tsv"
                echo "$block,$pos,$mode,$variant,$(cat "$tmp/seconds")"
            done
        done
    done
' bash "$tmp" "$fasta" "$oracle" > "$out/paired.csv"
{
    hostname -s
    gcc --version | head -1
    lscpu | awk '/^Model name:/'
    echo "threads=$OMP_NUM_THREADS"
    (cd "$tmp" && sha256sum previous final)
    sha256sum cmake/patches/wfa2-portability.cmake bench/cpu/alternatives/wfa_trace_{batch.cpp,ab.sbatch}
    sha256sum "$fasta" | awk '{print $1,"cohort.fa"}'
    sha256sum "$oracle" | awk '{print $1,"matrix.tsv"}'
} > "$out/environment.txt"
