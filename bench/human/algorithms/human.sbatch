#!/bin/bash -l
#SBATCH --job-name=human-biwfa
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
tmp=$(mktemp -d "${TMPDIR:-/tmp}/human-biwfa.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
build="$tmp/build"
out="results/human/algorithms/$SLURM_JOB_ID"
here="bench/human/algorithms"
mkdir -p "$out"
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-march=native \
      -DCMAKE_CXX_FLAGS=-march=native -DBUILD_TESTING=ON > "$tmp/config.log"
cmake --build "$build" --parallel 4 > "$tmp/build.log"
ctest --test-dir "$build" --output-on-failure > "$out/tests.txt"
flags=(-O3 -DNDEBUG -march=native -std=c++20 -fopenmp -Iinclude)
libs=("$build/libdistances.a" "$build/libedlib_static.a" "$build/libwfa_static.a" -lm)
g++ "${flags[@]}" "$here/trial.cpp" src/io/fasta.cpp "${libs[@]}" -o "$tmp/trial"
gcc -O3 -march=native -I"$build/_deps/wfa2-src" -c "$here/dispatch.c" -o "$tmp/dispatch.o"
g++ "${flags[@]}" "$here/trial.cpp" src/io/fasta.cpp "$tmp/dispatch.o" \
    -Wl,--wrap=wavefront_bialign,--wrap=wavefront_bialigner_set_sequences_bounds \
    "${libs[@]}" -o "$tmp/probe"
sha256sum CMakeLists.txt cmake/patches/wfa2-portability.cmake src/libraries/baselines.cpp \
    include/chps/baselines.hpp src/cli/main.cpp src/libraries/wfa2/wfa_adapter.c \
    tests/distances.cpp "$here"/{dispatch.c,study.py,trial.cpp,human.sbatch} > "$out/sources.sha256"
(cd "$tmp" && sha256sum trial probe build/editdist) > "$out/binaries.sha256"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
{ hostname -s; lscpu | grep -E 'Model name:|Thread.s. per core:';
  gcc --version | head -1; python3 -V; echo 'C/C++ flags: -O3 -DNDEBUG -march=native'; } > "$out/environment.txt"
srun --mpi=none --exact -N1 -n1 -c1 --cpu-bind=cores bash -c '
    sed -n "/^Cpus_allowed_list/p" /proc/self/status >> "$2/environment.txt"
    exec python3 bench/human/algorithms/study.py "$1" "$2" "$3"
' bash "$build" "$out" "$tmp"
