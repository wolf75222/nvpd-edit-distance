#!/bin/bash -l
#SBATCH --job-name=chps-oneview
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=maqao-%j.log
set -euo pipefail
phase=${1:?supply the folder containing src/cpu/dp.cpp}
inputs=${2:?supply prepared pair FASTAs}
out=${3:?supply a new output directory}
mkdir "$out"
out=$(realpath "$out")
inputs=$(realpath "$inputs")
source_root=$(realpath "$phase")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-oneview.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
maqao=${MAQAO:-maqao}
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
flags=(-O3 -g -march=native -fopenmp -std=c++20 -Wall -Wextra -Werror)
g++ "${flags[@]}" -I"$source_root/include" "$source_root/src/cpu/dp.cpp" "$source_root/src/io/fasta.cpp" \
    "$PWD/bench/cpu/maqao/kernel.cpp" -o "$work/kernel"
{
    hostname -s
    g++ --version | head -1
    "$maqao" --version | head -1
    echo "flags=${flags[*]}"
    echo "cpu_model=$(awk -F: '/model name/{print $2;exit}' /proc/cpuinfo)"
    sha256sum "$work/kernel" "$source_root/src/cpu/dp.cpp" "$source_root/src/cpu/dp_internal.hpp" bench/cpu/maqao/kernel.cpp
    if read -r paranoid < /proc/sys/kernel/perf_event_paranoid; then echo "perf_event_paranoid=$paranoid"; fi
    srun --cpu-bind=cores bash -c 'grep Cpus_allowed_list /proc/self/status'
} > "$out/environment.txt"
for threads in 1 4; do
    export OMP_NUM_THREADS=$threads
    xp="$work/threads-$threads"
    srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html,txt,csv \
        --lprof-params='sampling-rate=high' -xp="$xp" -- "$work/kernel" \
        "$inputs/near_32768.fa" 32768 "$threads" 16 > "$out/oneview-$threads.log" 2>&1
    tar -C "$xp" -czf "$out/threads-$threads.tar.gz" .
    echo "profile ready: $threads threads"
done
