#!/bin/bash -l
#SBATCH --job-name=nvpd-grace-placement
#SBATCH --partition=instant
#SBATCH --constraint=armgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
set -u
out=$(realpath -m "${1:?output directory}")
previous=$(realpath "${2:?completed scaling output}")
if (( ${SLURM_RESTART_COUNT:-0} > 0 )); then
    out="$out/resume-$SLURM_RESTART_COUNT"
fi
mkdir "$out"
cp "$previous/inputs.tsv" "$out/inputs.tsv"
work=$(mktemp -d "${TMPDIR:-/dev/shm}/chps-placement.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
bench="$PWD/bench/repeated/divergence/launch"
cmake -S . -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCHPS_WARNINGS_AS_ERRORS=ON -DCMAKE_PROJECT_INCLUDE="$bench/placement.cmake" >"$out/configure.log"
cmake --build "$work/build" --parallel 4 --target grace_placement >"$out/build.log"
{
    hostname -s
    printf 'job=%s\n' "$SLURM_JOB_ID"
    sha256sum "$bench"/placement.{cpp,cmake,sbatch,py} "$work/build/grace_placement"
    cat "$previous/inputs.sha256"
} >"$out/environment.txt"
srun --mpi=none --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
    python3 "$bench/placement.py" "$work/build/grace_placement" "$previous" "$out"
