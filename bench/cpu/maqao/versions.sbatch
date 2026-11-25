#!/bin/bash -l
#SBATCH --job-name=chps-maqao-versions
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=maqao-versions-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
if [[ -n ${MAQAO_HWLOC:-} ]]; then spack load "$MAQAO_HWLOC"; fi
set -u
original=$(realpath "${1:?source tree before the permutation}")
retained=$(realpath "${2:?source tree with the permutation}")
pair=$(realpath "${3:?near_32768.fa}")
mkdir "${4:?new output directory}"
out=$(realpath "$4")
rate=${5:-high}
repeats=${MAQAO_REPEATS:-96}
versions=(2025.1.4 2026.1.1)
if (( $# > 5 )); then versions=("${@:6}"); fi
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-maqao-versions.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close OMP_WAIT_POLICY=ACTIVE
driver="$PWD/bench/cpu/maqao/kernel.cpp"
flags=(-O3 -g -march=native -std=c++20 -fopenmp -Wall -Wextra -Werror)
if [[ ${MAQAO_FRAME_POINTERS:-0} == 1 ]]; then flags+=(-fno-omit-frame-pointer); fi
for backend in original retained; do
    tree=$original
    [[ $backend == original ]] || tree=$retained
    g++ "${flags[@]}" -I"$tree/include" "$tree/src/cpu/dp.cpp" "$tree/src/io/fasta.cpp" "$driver" -o "$work/$backend"
    "$work/$backend" "$pair" 32768 "$OMP_NUM_THREADS" 1 | awk -F, 'NR == 2 && $1 == 1479 {ok = 1} END {exit !ok}'
done
cd "$work"
hostname -s > "$out/environment.txt"
g++ --version | sed -n '1p' >> "$out/environment.txt"
{
    printf 'job=%s\nthreads=%s\nsampling_rate=%s\nrepeats=%s\n' "$SLURM_JOB_ID" "$OMP_NUM_THREADS" "$rate" "$repeats"
    printf 'flags=%s\n' "${flags[*]}"
    for param in perf_event_mlock_kb perf_event_paranoid; do
        if [[ -r /proc/sys/kernel/$param ]]; then printf '%s=%s\n' "$param" "$(cat "/proc/sys/kernel/$param")"; fi
    done
    command -v lstopo-no-graphics || true
} >> "$out/environment.txt"
sha256sum original retained >> "$out/environment.txt"
for version in "${versions[@]}"; do
    maqao="$(spack location -i "maqao@$version")/bin/maqao"
    "$maqao" --version | sed -n '1p' >> "$out/environment.txt"
    for backend in original retained; do
        xp="$work/$version-$backend"
        srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html \
            --lprof-params="sampling-rate=$rate" -xp="$xp" -- "$work/$backend" "$pair" 32768 "$OMP_NUM_THREADS" "$repeats" \
            > "$out/$version-$backend.log" 2>&1
        tar -czf "$out/$version-$backend.tar.gz" -C "$xp" .
    done
done
