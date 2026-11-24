#!/bin/bash -l
#SBATCH --job-name=chps-chr21
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=16G
#SBATCH --constraint=armgpu
#SBATCH --time=00:45:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE
build=$(mktemp -d "${TMPDIR:-/tmp}/chr21.XXXXXXXX")
trap cleanup EXIT
cleanup() { rm -rf -- "$build"; }
out=results/human/$SLURM_JOB_ID
mkdir -p "$out"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv > "$out/device.txt"
nvcc --version >> "$out/device.txt"
lscpu | awk -F: '/^Architecture:|^Model name:/{print}' >> "$out/device.txt"
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 -DBUILD_TESTING=ON
cmake --build "$build" --parallel "$SLURM_CPUS_PER_TASK"
step=(srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores)
"${step[@]}" ctest --test-dir "$build" --output-on-failure | sed "s|$build|build|g" > "$out/validation.txt"
"${step[@]}" python3 bench/human_run.py --binary "$build/editdist" --threads "$SLURM_CPUS_PER_TASK" --output "$out/trials.csv"
if [[ -n ${KOKKOS_SOURCE:-} ]]; then
    export CHPS_HUMAN_PAIR="$PWD/data/human/chr21-2m.fa"
    CHPS_HUMAN_EXPECTED=$(python3 -c 'import csv,sys; print(next(r["distance"] for r in csv.DictReader(open(sys.argv[1])) if r["case"]=="chr21-2m" and r["algorithm"]=="edlib" and r["status"]=="ok"))' "$out/trials.csv")
    export CHPS_HUMAN_EXPECTED
    (
        set +u
        spack unload /dvle736
        spack load /uxlsycg
        set -u
        "${step[@]}" timeout 420 bash bench/kokkos-step.sh gpu "$out/kokkos"
    )
fi
