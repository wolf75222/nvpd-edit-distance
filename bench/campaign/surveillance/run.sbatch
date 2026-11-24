#!/bin/bash -l
#SBATCH --job-name=nvpd-cohort
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=cohort-%j.log
set -eo pipefail
target=$1
cuda=(-DCHPS_CUDA=OFF)
case "$target" in
    cpu|release) romeo_load_x64cpu_env >/dev/null; spack load /q22mihh ;;
    gpu) romeo_load_armgpu_env >/dev/null; spack load /dvle736; spack load /phhwuh3
         cuda=(-DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90) ;;
    *) exit 2 ;;
esac
set -u
out="results/phylogeny/sars-romeo/$SLURM_JOB_ID"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-cohort.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$out"
cmake -S . -B "$tmp/app" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native "${cuda[@]}" > "$tmp/configure.log"
cmake --build "$tmp/app" --parallel "$SLURM_CPUS_PER_TASK" --target editdist
{
    hostname -s
    lscpu | grep -E 'Model name|Architecture:|Socket\(s\):|Core\(s\) per socket:'
    g++ --version | head -1
    if [[ $target == gpu ]]; then nvidia-smi --query-gpu=name,driver_version --format=csv; fi
    sha256sum "$2" "$3" | cut -d' ' -f1
} > "$out/environment.txt"
(cd "$tmp/app" && sha256sum editdist) > "$out/binary.sha256"
sha256sum include/chps/*.hpp src/{cli/main.cpp,cpu/*.cpp,cuda/*.cu} bench/campaign/surveillance/{measure.py,run.sbatch} > "$out/sources.sha256"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false
srun --ntasks=1 --cpu-bind=cores --cpus-per-task="$SLURM_CPUS_PER_TASK" python3 \
    bench/campaign/surveillance/measure.py "$tmp/app/editdist" "$2" "$3" "$out/trials.csv" "$target"
