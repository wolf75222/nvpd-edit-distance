#!/bin/bash -l
#SBATCH --job-name=nvpd-e2e
#SBATCH --partition=instant
#SBATCH --constraint=armgpu
#SBATCH --nodelist=romeo-a002
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-socket=1
#SBATCH --cores-per-socket=1
#SBATCH --cpus-per-task=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem=8G
#SBATCH --time=00:58:00
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
set -u
root=${CHPS_ROOT:-$PWD}
out=$(realpath -m "${1:?output directory}")
inputs=$(realpath "${2:?four-chromosome manifest}")
study="$root/bench/human/grace-numa/amortization"
mkdir "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/grace-amortization.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
{
    hostname -s
    date -Iseconds
    scontrol show job -dd "$SLURM_JOB_ID" | awk '/JobId=|CPU_IDs=|NumCPUs=|AllocTRES=|Partition=/'
    g++ --version | head -1
    sha256sum "$study"/{amortization.cpp,amortization.cmake,measure_e2e.py,cache_probe.py,run.sbatch}
} >"$out/environment.txt"
cmake -S "$root" -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCHPS_WARNINGS_AS_ERRORS=ON -DCMAKE_PROJECT_INCLUDE="$study/amortization.cmake" >"$out/configure.log"
cmake --build "$work/build" --parallel 2 --target grace_amortization >"$out/build.log"
sha256sum "$work/build/grace_amortization" >>"$out/environment.txt"
taskset -c 0,72 python3 "$study/cache_probe.py" "$inputs" "$out/cache.csv" initial
panel=$(dirname "$(head -1 "$inputs" | cut -f1)")
(cd "$panel" && sha256sum --check "$root/results/human/grace-scaling/inputs.sha256")
# Match every manifest path to the file just checked.
python3 - "$inputs" "$panel" <<'PY'
import sys
from pathlib import Path
manifest, panel = map(Path, sys.argv[1:])
paths = [Path(row.split("\t")[0]).resolve() for row in manifest.read_text().splitlines()]
if paths != [(panel / f"chr{i}.fa").resolve() for i in range(1, 5)]:
    raise ValueError("Expected the four checked FASTA files from CHM13")
PY
cp "$root/results/human/grace-scaling/inputs.sha256" "$out/"
taskset -c 0,72 python3 "$study/cache_probe.py" "$inputs" "$out/cache.csv" after_hash
taskset -c 0,72 python3 -u "$study/measure_e2e.py" "$work/build/grace_amortization" "$inputs" "$out"
taskset -c 0,72 python3 "$study/cache_probe.py" "$inputs" "$out/cache.csv" after_trials
