#!/bin/bash -l
#SBATCH --job-name=nvpd-repeat
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --output=repeat-%j.log
set -eo pipefail
platform=${1:?expected epyc, epyc-trace, epyc-bands, epyc-scalar or gh200}
if [[ $platform == epyc* ]]; then
    romeo_load_x64cpu_env >/dev/null
    spack load /q22mihh
elif [[ $platform == gh200 ]]; then
    romeo_load_armgpu_env >/dev/null
    spack load /phhwuh3
    spack load /dvle736
else
    exit 2
fi
set -u
build=$(mktemp -d "${TMPDIR:-/tmp}/repeat-nvpd.XXXXXX")
trap 'rm -rf -- "$build"' EXIT
out=results/repeated/jobs/$SLURM_JOB_ID
mkdir -p "$out"
opts=(-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON "-DCMAKE_PROJECT_INCLUDE=$PWD/bench/repeated/targets.cmake")
if [[ $platform == epyc ]]; then
    if [[ -n ${CHPS_REPEAT_INITIAL:-} ]]; then cp "$CHPS_REPEAT_INITIAL" "$build/initial.cpp"
    else git show f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83:src/dp.cpp > "$build/initial.cpp"; fi
    echo "95a75e9e3dc1d1540aafa73ba05313269fdfbcd3740453d169a3c041446d63c6  $build/initial.cpp" | sha256sum -c --status
    opts+=("-DCHPS_REPEAT_INITIAL=$build/initial.cpp" -DCMAKE_CXX_FLAGS=-march=native)
elif [[ $platform == epyc-* ]]; then
    opts+=(-DCMAKE_CXX_FLAGS=-march=native)
else
    opts+=(-DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90)
fi
cmake -S . -B "$build" "${opts[@]}"
cmake --build "$build" --parallel "$(( SLURM_CPUS_PER_TASK < 4 ? SLURM_CPUS_PER_TASK : 4 ))"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
if [[ $platform == epyc ]]; then export OMP_NUM_THREADS=8; else export OMP_NUM_THREADS=1; fi
ctest --test-dir "$build" --output-on-failure | sed "s|$build|build|g" > "$out/tests.txt"
{
  hostname -s
  g++ --version | head -1
  cmake --version | head -1
  lscpu | awk '/^(Architecture|Model name):/'
  if [[ $platform == gh200 ]]; then nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv; nvcc --version | tail -2; fi
  env | sort | sed -n '/^OMP_/p'
  cat "$build/CMakeFiles/distances.dir/flags.make" | sed -n '/^CXX_FLAGS =/p;/^CUDA_FLAGS =/p'
  sha256sum src/cpu/dp.cpp src/cpu/dp_internal.hpp src/cpu/banded.cpp src/libraries/baselines.cpp src/cuda/cuda.cu src/cpu/hirschberg.cpp bench/repeated/{inputs.csv,kernel.cpp,measure.py,run.sbatch,targets.cmake}
} > "$out/environment.txt"
export build out platform
if [[ $platform == epyc-bands ]]; then
    sha256sum bench/repeated/bands.py >> "$out/environment.txt"
    srun --mpi=none --cpu-bind=cores python3 bench/repeated/bands.py "$build/repeat_kernel" "${BAND_INPUTS:-inputs}" "$out"
    exit 0
fi
srun --mpi=none --cpu-bind=cores python3 bench/repeated/divergence/measure.py "$build" "$out/divergence.csv"
