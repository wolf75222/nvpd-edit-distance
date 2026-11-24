#!/bin/bash -l
#SBATCH --job-name=chps-read-blocks
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:25:00
#SBATCH --output=reader-%j.log
set -eu -o pipefail
data=$(realpath "${1:-data/human/chm13-panel}")
out=${PANEL_OUTPUT:-results/human/chromosome-panel/$SLURM_JOB_ID}
mkdir -p "$out"
out=$(realpath "$out")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-reader.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
if [[ -n ${PANEL_PORTABLE:-} ]]; then
    cp "$PANEL_PORTABLE" "$work/portable.cpp"
else
    git show 421927e9:src/io/fasta.cpp > "$work/portable.cpp"
fi
echo "2bb6087a2ae06544a0c20185eafa06b53b7b755a1a5a41744f66d9fa2b685539  $work/portable.cpp" | sha256sum -c --status
bench=bench/human/chromosome-panel/blocks
for variant in portable blocks; do
    source="$work/portable.cpp"
    if [[ $variant == blocks ]]; then source=src/io/fasta.cpp; fi
    g++ -std=c++20 -O3 -march=native -fopenmp -Iinclude "$source" "$bench/probe.cpp" -o "$work/$variant"
done
{
    hostname -s
    g++ --version | head -1
    sha256sum src/io/fasta.cpp "$bench"/{probe.cpp,sweep.py,run.sbatch} "$bench/../counters.hpp"
} > "$out/environment.txt"
srun --cpu-bind=cores --distribution=block:block python3 "$bench/sweep.py" "$data" "$work" "$out"
