#!/bin/bash -l
#SBATCH --job-name=nvpd-adaptive
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:55:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-adaptive.XXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
build="$tmp/build"
out="results/adaptive/$SLURM_JOB_ID"
mkdir -p "$out"
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-march=native \
      -DCMAKE_CXX_FLAGS=-march=native -DBUILD_TESTING=ON > "$tmp/config.log"
cmake --build "$build" --parallel 4 > "$tmp/build.log"
ctest --test-dir "$build" --output-on-failure > "$out/tests.txt"
g++ -O3 -DNDEBUG -march=native -std=c++20 -fopenmp -Iinclude -Isrc/cpu \
    bench/cpu/adaptive/driver.cpp src/io/fasta.cpp "$build/libdistances.a" "$build/libedlib_static.a" \
    "$build/libwfa_static.a" -lm -o "$tmp/driver"
g++ -O3 -DNDEBUG -march=native -std=c++20 -fopenmp -Iinclude -Isrc/cpu \
    -fopt-info-vec-optimized="$tmp/vectorization" -c src/cpu/nvpd_adaptive.cpp -o "$tmp/vector.o"
# The control differs only in the C++ spelling of a three-way minimum. GCC
# 11.5 scalarizes the indirect load with its initializer-list overload.
mkdir "$tmp/control"
cp src/cpu/nvpd_adaptive.cpp "$tmp/control/"
sed 's/std::min(inf, std::min(1 + std::min(prev\[p\], prev\[p+1\]), via_match))/std::min({inf, 1 + std::min(prev[p], prev[p+1]), via_match})/' \
    src/cpu/nvpd_band_row.hpp > "$tmp/control/nvpd_band_row.hpp"
cmp -s src/cpu/nvpd_band_row.hpp "$tmp/control/nvpd_band_row.hpp" && exit 1
g++ -O3 -DNDEBUG -march=native -std=c++20 -fopenmp -Iinclude -Isrc/cpu \
    -Dnvpd_band_row=nvpd_band_row_control -fopt-info-vec-all="$tmp/control/report" -c "$tmp/control/nvpd_adaptive.cpp" -o "$tmp/control.o"
g++ -O3 -DNDEBUG -march=native -std=c++20 -fopenmp -Iinclude -Isrc/cpu \
    bench/cpu/adaptive/driver.cpp src/io/fasta.cpp "$tmp/control.o" "$build/libdistances.a" \
    "$build/libedlib_static.a" "$build/libwfa_static.a" -lm -o "$tmp/driver-control"
{ grep 'nvpd_band_row.hpp' "$tmp/vectorization" | sort -u
  objdump -d "$tmp/vector.o" | awk '/vpgatherdd/ && ++n<=2 {print "SIMD gather", $0}'
  grep -E 'not suitable for gather|couldn.t vectorize loop' "$tmp/control/report" | \
      grep nvpd_band_row.hpp | sed -E 's|.*/nvpd_band_row.hpp|scalar control nvpd_band_row.hpp|' | sed -n '1,2p'; } > "$out/vectorization.txt"
(cd "$tmp" && sha256sum driver driver-control control/nvpd_band_row.hpp) > "$out/binaries.sha256"
sha256sum CMakeLists.txt src/cpu/{dp.cpp,dp_internal.hpp,nvpd_adaptive.cpp,nvpd_band_row.hpp} \
    src/libraries/{baselines.cpp,wfa2/wfa_adapter.c} include/chps/adaptive.hpp \
    src/io/fasta.cpp bench/cpu/adaptive/{driver.cpp,cases.py,measure.py,run.sbatch} > "$out/sources.sha256"
perf="${PERF_BIN:-/project/$SLURM_JOB_ACCOUNT/$USER/tools/perf/bin/perf}"
{ hostname -s; lscpu | grep -E 'Model name:|Thread.s. per core:'; gcc --version | head -1;
  python3 -V; "$perf" --version;
  if IFS= read -r value < /proc/sys/kernel/perf_event_paranoid; then echo "perf_event_paranoid=$value"; fi; } > "$out/environment.txt"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
srun --mpi=none --exact -N1 -n1 -c1 --cpu-bind=cores bash -c '
    sed -n "/^Cpus_allowed_list/p" /proc/self/status >> "$3/environment.txt"
    for c in /sys/devices/system/cpu/cpu0/cache/index*; do paste -sd " " "$c"/{level,type,size,coherency_line_size} | xargs; done >> "$3/environment.txt"
    exec python3 bench/cpu/adaptive/measure.py "$@"
' bash "$build" "$tmp" "$out" "$perf"
