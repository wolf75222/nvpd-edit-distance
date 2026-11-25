#!/bin/bash -l
#SBATCH --job-name=chps-hg002-flanks
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=10G
#SBATCH --time=00:14:00
#SBATCH --output=hg002-flanks-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
inputs=$(realpath "${1:?HG002 chromosome FASTA directory}")
out=$(realpath -m "${2:?new output directory}")
mkdir "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-flanks.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
if [[ -n ${3:-} ]]; then
    cp "$3" "$work/split-before.cpp"
else
    git show 0c916d5:src/cpu/split.cpp > "$work/split-before.cpp"
fi
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close OMP_WAIT_POLICY=PASSIVE
cmake -S . -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG -march=native" \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native" -DBUILD_TESTING=ON -DCHPS_WARNINGS_AS_ERRORS=ON \
    -DCMAKE_PROJECT_INCLUDE="$PWD/bench/cpu/chromosome-dispatch/flanks/targets.cmake" \
    -DCHPS_SPLIT_BASELINE="$work/split-before.cpp"
cmake --build "$work/build" --target flank_trial check_split -j 8
ctest --test-dir "$work/build" -R '^split$' --output-on-failure
{
    hostname -s
    g++ --version | head -1
    lscpu | awk '/^Model name:/'
    srun --mpi=none --cpu-bind=cores bash -c 'grep Cpus_allowed_list: /proc/self/status'
    sha256sum "$work/build/flank_trial" src/cpu/split.cpp "$work/split-before.cpp"
} > "$out/environment.txt"
python3 - "$inputs" <<'PY'
from pathlib import Path
import csv,hashlib,sys
root=Path(sys.argv[1])
with open('results/human/chromosome-dispatch/manifest.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        path=root/f"{r['chrom']}.v{r['version']}.fa"
        h=hashlib.sha256()
        with path.open('rb') as stream:
            for line in stream:
                if not line.startswith(b'>'): h.update(line.rstrip())
        if h.hexdigest()!=r['sha256_bases']:
            raise ValueError(f'Unexpected bases in {path}')
PY
while IFS=, read -r chrom _ _ score _; do
    [[ $chrom != chrom ]] || continue
    srun --input=none --mpi=none --cpu-bind=cores "$work/build/flank_trial" "$chrom" \
       "$inputs/$chrom.v1.1.fa" "$inputs/$chrom.v1.2.fa" "$score" > "$out/$chrom.csv"
done < results/human/chromosome-dispatch/partitions.csv

python3 - "$out" <<'PYTHON'
from pathlib import Path
import csv, sys
out = Path(sys.argv[1])
for chrom in ("chr1_MATERNAL","chr1_PATERNAL","chr21_MATERNAL","chr21_PATERNAL"):
    with (out/f"{chrom}.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 64 or len({(r["block"],r["method"],r["threads"]) for r in rows}) != 64:
        raise ValueError(f"Incomplete timing series for {chrom}")
PYTHON
