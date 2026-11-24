#!/usr/bin/env bash
# Run inside one allocated CPU/GPU task in the project root.
set -euo pipefail
target=${1:?specify cpu or gpu}
out=${2:?specify output directory}
workers=${KOKKOS_THREADS:-${SLURM_CPUS_PER_TASK:-4}}
(( workers <= 8 )) || workers=8
mkdir -p "$out"
out=$(realpath "$out")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/chps-kokkos.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
pin=175257a51ff29a0059ec48bcd233ee096b2c0438
host_compiler=${CXX:-c++}
if [[ -n ${KOKKOS_PREFIX:-} ]]; then
    prefix=$KOKKOS_PREFIX
    libversion=external
    compiler=${CXX:-c++}
else
    : "${KOKKOS_SOURCE:?set to a Kokkos Git checkout containing tag 4.5.01}"
    src=$tmp/src
    mkdir "$src"
    git -C "$KOKKOS_SOURCE" archive "$pin" | tar -xC "$src"
    prefix=$tmp/install
    libversion="4.5.01 $pin"
    options=(-DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_TESTS=OFF -DKokkos_ENABLE_EXAMPLES=OFF)
    compiler=${CXX:-c++}
    if [[ $target == gpu ]]; then
        compiler=$src/bin/nvcc_wrapper
        # The pinned wrapper targets sm_70 before CMake selects an architecture.
        # Configure Hopper in the early compiler probes as well as the kernels.
        sed -i 's/^default_arch="sm_70"$/default_arch="sm_90"/' "$compiler"
        options+=(-DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON)
    fi
    cmake -S "$src" -B "$tmp/klib" -DCMAKE_CXX_STANDARD=20 -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_CXX_COMPILER="$compiler" "${options[@]}"
    cmake --build "$tmp/klib" --target install --parallel "$workers"
fi
opts=(-DCHPS_KOKKOS=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$prefix")
if [[ $target == gpu ]]; then
    opts+=(-DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90)
elif [[ $target == cpu ]]; then
    opts+=(-DCHPS_CUDA=OFF "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -march=native")
else
    echo 'target must be cpu or gpu' >&2
    exit 2
fi
# Kokkos's CMake launcher selects nvcc_wrapper for dependent targets only.
# Third-party host code such as edlib must retain its ordinary C++ compiler.
cmake -S . -B "$tmp/app" -DCMAKE_CXX_COMPILER="$host_compiler" "${opts[@]}"
cmake --build "$tmp/app" --target editdist measure_kokkos check_distances --parallel "$workers"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_NUM_THREADS="$workers"
ctest --test-dir "$tmp/app" --output-on-failure
{
    echo "Kokkos $libversion"
    echo "target=$target host_workers=$workers"
    echo "job=${SLURM_JOB_ID:-local}"
    "$compiler" --version | head -2
    uname -m
} > "$out/environment.txt"
python3 bench/kokkos/kokkos_measure.py "$target" "$tmp/app" "$tmp" "$out" "$workers"
