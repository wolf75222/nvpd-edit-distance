#!/bin/bash
set -euo pipefail
target=$1
shift
build=$(mktemp -d "${TMPDIR:-/tmp}/chps-edit-${SLURM_JOB_ID}.XXXXXXXX")
trap 'rm -rf -- "$build"' EXIT
mkdir -p results
cuda_args=(-DCHPS_CUDA=OFF)
if [[ $target == gpu ]]; then
    cuda_args=(-DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90)
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
    nvcc --version
fi
hostname
lscpu | awk -F: '/^Model name:|^Architecture:/{print}'
cmake --version | head -n 1
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON "${cuda_args[@]}"
cmake --build "$build" --parallel "$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE
step=(srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores)
"${step[@]}" ctest --test-dir "$build" --output-on-failure
status=0
if [[ -n ${COHORT_FASTA:-} ]]; then
    groups=()
    if [[ -n ${COHORT_GROUPS:-} ]]; then groups=(--groups "$COHORT_GROUPS"); fi
    if [[ -n ${COHORT_VERIFY_WITH:-} ]]; then groups+=(--verify-with "$COHORT_VERIFY_WITH"); fi
    "${step[@]}" python3 bench/campaign/cohort.py --binary "$build/editdist" --fasta "$COHORT_FASTA" \
        --threads "$SLURM_CPUS_PER_TASK" --output "results/cohort-${SLURM_JOB_ID}" "${groups[@]}" || status=$?
fi
"${step[@]}" python3 bench/campaign/run.py --target "$target" --binary "$build/editdist" \
    --output "results/${target}-${SLURM_JOB_ID}.csv" "$@" || status=$?
exit "$status"
