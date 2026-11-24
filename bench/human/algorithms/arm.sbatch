#!/bin/bash -l
#SBATCH --job-name=chr21q-biwfa
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
tmp=$(mktemp -d "${TMPDIR:-/tmp}/chr21q-biwfa.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
build="$tmp/build"
out="results/human/algorithms/arm/$SLURM_JOB_ID"
here="bench/human/algorithms"
mkdir -p "$out"
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-march=native \
      -DCMAKE_CXX_FLAGS=-march=native -DBUILD_TESTING=ON > "$tmp/config.log"
cmake --build "$build" --parallel 4 --target editdist check_distances > "$tmp/build.log"
ctest --test-dir "$build" --output-on-failure > "$out/tests.txt"
g++ -O3 -DNDEBUG -march=native -std=c++20 -fopenmp -Iinclude \
    "$here/trial.cpp" src/io/fasta.cpp "$build/libdistances.a" "$build/libedlib_static.a" \
    "$build/libwfa_static.a" -lm -o "$tmp/trial"
sha256sum CMakeLists.txt src/io/fasta.cpp cmake/patches/wfa2-portability.cmake src/libraries/baselines.cpp \
    include/chps/baselines.hpp src/cli/main.cpp src/libraries/wfa2/wfa_adapter.c \
    "$here"/{arm.py,arm.sbatch,study.py,trial.cpp} > "$out/sources.sha256"
(cd "$tmp" && sha256sum trial build/editdist) > "$out/binaries.sha256"
{ hostname -s; lscpu | grep -E 'Model name:|Thread.s. per core:';
  gcc --version | head -1; python3 -V; echo 'C/C++ flags: -O3 -DNDEBUG -march=native'; } > "$out/environment.txt"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1
srun --mpi=none --exact -N1 -n1 -c4 --cpu-bind=cores bash -c '
    sed -n "/^Cpus_allowed_list/p" /proc/self/status >> "$2/environment.txt"
    exec python3 bench/human/algorithms/arm.py "$@"
' bash "$build" "$out" "$tmp" "${CHPS_CHAIN:-data/human/hg38ToHs1.over.chain.gz}"
