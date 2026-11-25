#!/bin/bash -l
#SBATCH --job-name=chps-frame-cost
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=maqao-frames-%j.log
set -euo pipefail
source_root=$(realpath "${1:?retained source tree}")
inputs=$(realpath "${2:?directory containing the chromosome pairs and inputs.csv}")
out=${3:?new output CSV}
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-frames.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
flags=(-O3 -g -march=native -std=c++20 -fopenmp -Wall -Wextra -Werror)
for variant in retained frames; do
    extra=()
    [[ $variant != frames ]] || extra+=(-fno-omit-frame-pointer)
    g++ "${flags[@]}" "${extra[@]}" -I"$source_root/include" "$source_root/src/cpu/dp.cpp" \
        "$source_root/src/io/fasta.cpp" bench/cpu/maqao/kernel.cpp -o "$work/$variant"
done
srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$out" \
    "$work/retained" "$work/frames" --lengths 32768 100000 --threads 1 4 --blocks 8
