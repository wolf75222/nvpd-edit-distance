#!/bin/bash -l
#SBATCH --job-name=chps-gpu-sectors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --constraint=armgpu
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
set -u
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
build=$(mktemp -d "${TMPDIR:-/tmp}/chps-sectors.XXXXXXXX")
trap 'rm -rf -- "$build"' EXIT
out="results/locality/$SLURM_JOB_ID"
mkdir -p "$out"
{
    date -u +%FT%TZ
    hostname -s
    nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv
    nvcc --version
    ncu --version
} > "$out/device.txt"
sha256sum bench/gpu_sectors.{cu,sbatch} scripts/synthetic.py src/{cuda.cu,dp.cpp,dp_internal.hpp,fasta.cpp} > "$out/source.sha256"
nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -Iinclude -Isrc \
    bench/gpu_sectors.cu src/cuda.cu src/fasta.cpp src/dp.cpp -o "$build/gpu_sectors"
python3 scripts/synthetic.py --out "$build/data" --lengths 5000 --rates 0.02 --pairs 1 --seed 2026
cp "$build"/data/pairs/*.fa "$build/standard.fa"
python3 - "$build" <<'PY'
import sys
from pathlib import Path
work = Path(sys.argv[1])
a, b = ["".join(r.splitlines()[1:]) for r in (work/"standard.fa").read_text().split(">")[1:]]
# The score-only production kernel places the longer string in its target columns.
if len(a) > len(b):
    a, b = b, a
# Introduce a long interval without A after a known A.
rare = b[:128] + "A" + "C"*160 + b[289:]
wide = (b * (300000 // len(b) + 1))[:300000]
for name, left, right in [("rare", a, rare), ("wide", a[:32], wide)]:
    (work/(name+".fa")).write_text(f">a\n{left}\n>b\n{right}\n")
PY
(cd "$build" && sha256sum standard.fa rare.fa wide.fa) > "$out/inputs.sha256"
probe_metrics='l1tex__t_sector_hit_rate.pct,lts__t_sector_hit_rate.pct'
echo "$probe_metrics" > "$out/requested-metrics.txt"
if ncu --profile-from-start off --kernel-name 'regex:row_tile' --launch-count 1 \
        --metrics "$probe_metrics" --csv --page raw --log-file "$build/probe.csv" \
        "$build/gpu_sectors" capture standard "$build/standard.fa" > "$build/ncu.log" 2>&1; then
    cp "$build/probe.csv" "$out/cache-probe.csv"
else
    if ! grep -h ERR_NVGPUCTRPERM "$build/probe.csv" "$build/ncu.log" > "$out/cache-permission.txt"; then
        cat "$build/probe.csv" "$build/ncu.log" >&2
        exit 1
    fi
fi
echo 'case,schedule,kernel,rows,columns,reach,sm_count,output_width,launches,distance,access,warp_requests,lane_words,distinct_words,sectors32' > "$out/sectors.csv"
for case in standard rare wide; do
    "$build/gpu_sectors" model "$case" "$build/$case.fa" >> "$out/sectors.csv"
done
python3 - "$out/sectors.csv" <<'PYTHON'
import csv, sys
with open(sys.argv[1]) as f:
    observed = {r["case"]: (r["kernel"], int(r["output_width"]))
                for r in csv.DictReader(f) if r["schedule"] == "adaptive"}
expected = {"standard": ("row_tile", 256), "rare": ("next_row", 256), "wide": ("row_tile", 1024)}
if observed != expected:
    raise ValueError(f"cases did not exercise their intended dispatch: {observed}")
PYTHON
