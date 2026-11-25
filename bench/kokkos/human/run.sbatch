#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --job-name=kokkos-human
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=kokkos-human-%j.log
set -eo pipefail
target=${1:?expected epyc or gh200}
if [[ $target == epyc ]]; then
    romeo_load_x64cpu_env >/dev/null
    spack load /q22mihh
else
    romeo_load_armgpu_env >/dev/null
    spack load /phhwuh3
    spack load /uxlsycg
fi
set -u
tmp=$(mktemp -d "${TMPDIR:-/tmp}/kokkos-human.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
out=$PWD/results/kokkos/human/jobs/$SLURM_JOB_ID
mkdir -p "$out" "$tmp/ksrc"
pin=175257a51ff29a0059ec48bcd233ee096b2c0438
git -C "${KOKKOS_SOURCE:?set path to pinned Kokkos checkout}" archive "$pin" | tar -xC "$tmp/ksrc"
install=$tmp/install
host=${CXX:-c++}
compiler=$host
ko=(-DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_TESTS=OFF -DKokkos_ENABLE_EXAMPLES=OFF)
app=(-DCHPS_KOKKOS=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$install"
     "-DCMAKE_PROJECT_INCLUDE=$PWD/bench/kokkos/human/targets.cmake" -DCMAKE_CXX_COMPILER="$host")
if [[ $target == gh200 ]]; then
    compiler=$tmp/ksrc/bin/nvcc_wrapper
    sed -i 's/^default_arch="sm_70"$/default_arch="sm_90"/' "$compiler"
    ko+=(-DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON)
    app+=(-DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90)
else
    app+=(-DCHPS_CUDA=OFF "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -march=native")
fi
if [[ ${KOKKOS_PROFILE:-0} == 1 ]]; then app+=(-DCHPS_HUMAN_PROFILE=ON); fi
cmake -S "$tmp/ksrc" -B "$tmp/kbuild" -DCMAKE_CXX_STANDARD=20 -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$install" -DCMAKE_CXX_COMPILER="$compiler" "${ko[@]}"
cmake --build "$tmp/kbuild" --target install --parallel 4
cmake -S . -B "$tmp/app" "${app[@]}"
cmake --build "$tmp/app" --parallel 4
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE OMP_NUM_THREADS=2
ctest --test-dir "$tmp/app" --output-on-failure | sed "s|$tmp|build|g" > "$out/validation.txt"
"$tmp/app/human_kokkos" --check
echo '101 differential pairs, cost + trace, reference and candidate' >> "$out/validation.txt"
{
  hostname -s
  lscpu | awk '/^(Architecture|Model name):/'
  echo "Kokkos 4.5.01 $pin"
  "$host" --version | head -1
  if [[ $target == gh200 ]]; then nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv; nvcc --version | tail -2; fi
  sed -n '/^CXX_FLAGS =/p;/^CUDA_FLAGS =/p' "$tmp/app/CMakeFiles/distances.dir/flags.make"
  sha256sum src/kokkos/{kokkos.cpp,row_policy.hpp} bench/kokkos/human/{inputs.py,trial.cpp,measure.py,targets.cmake,run.sbatch}
} > "$out/environment.txt"
inputs=${KOKKOS_INPUTS:-inputs}
python3 bench/kokkos/human/measure.py reference "$tmp/app/editdist" "$inputs" "$out"
if [[ ${KOKKOS_PROFILE:-0} == 1 ]]; then
    nsys --version >> "$out/environment.txt"
    srun --mpi=none --cpu-bind=cores python3 bench/kokkos/human/profile.py "$tmp/app/human_kokkos" "$inputs" "$out" "$tmp"
    exit 0
fi
srun --mpi=none --cpu-bind=cores --distribution=block:block python3 bench/kokkos/human/measure.py \
    "$target" "$tmp/app/human_kokkos" "$inputs" "$out"
