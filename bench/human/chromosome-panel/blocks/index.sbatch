#!/bin/bash -l
#SBATCH --job-name=chps-fai
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=reader-index-%j.log
set -euo pipefail
data=$(realpath "${1:-data/human/chm13-panel}")
out=${PANEL_OUTPUT:-results/human/chromosome-panel/index-$SLURM_JOB_ID}
mkdir -p "$out"
out=$(realpath "$out")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-fai.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
base=bench/human/chromosome-panel
flags=(-std=c++20 -O3 -march=native -fopenmp -DCHPS_PANEL_FAI -Iinclude)
g++ "${flags[@]}" src/io/fasta.cpp "$base/blocks/probe.cpp" -o "$work/probe"
{
    hostname -s
    lscpu | awk '/^(Architecture|Model name|Thread\(s\) per core):/'
    g++ --version | head -1
    echo "g++ ${flags[*]}"
    sha256sum src/io/fasta.cpp include/chps/fasta.hpp "$base/counters.hpp" "$base/measure.py" "$base/blocks/"{probe.cpp,index.py,index.sbatch}
    sha256sum "$work/probe" | awk '{print $1, "probe"}'
} > "$out/environment.txt"
srun --exclusive --cpu-bind=cores --distribution=block:block python3 "$base/blocks/index.py" "$data" "$work/probe" "$out"
