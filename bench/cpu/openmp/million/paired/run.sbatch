#!/bin/bash -l
#SBATCH --job-name=chps-paired-million
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=6G
#SBATCH --time=00:35:00
#SBATCH --output=million-eight-%j.log
set -euo pipefail
data=$(realpath "${1:?synthetic FASTA root}")
out=$(realpath -m "${2:?new output directory}")
mkdir "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-eight.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
g++ -std=c++20 -fopenmp -O3 -g -march=native -Wall -Wextra -Werror -Iinclude \
    bench/cpu/openmp/million/paired/trial.cpp src/io/fasta.cpp src/cpu/dp.cpp -o "$work/trial"
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=8
export out data work
{
    hostname -s
    g++ --version | head -1
    lscpu | awk '/^Model name:|Thread\(s\) per core:|Core\(s\) per socket:/'
    sha256sum "$work/trial" src/cpu/{dp.cpp,nvpd_avx512.hpp}
} > "$out/environment.txt"
python3 - "$data" "$out" <<'PY'
from pathlib import Path
import csv, hashlib, sys
data, out = map(Path, sys.argv[1:])
with open('results/romeo/completion-1m/references.csv') as f:
    rows = [r for r in csv.DictReader(f) if r['rate']=='0.05' and r['length']=='1000000']
with (out/'manifest.tsv').open('w') as f:
    print('pair','path','sha256','score',sep='\t',file=f)
    for i in range(1,6):
        refs=[r for r in rows if int(r['rep'])==i]
        r=refs[0]; path=data/r['path']
        if ({x['algorithm'] for x in refs}!={'edlib','wfa2'} or len({(x['distance'],x['path'],x['input_sha256']) for x in refs})!=1
                or hashlib.sha256(path.read_bytes()).hexdigest()!=r['input_sha256']):
            raise ValueError(f'Oracle or input mismatch for pair {i}')
        print(i,r['path'],r['input_sha256'],r['distance'],sep='\t',file=f)
PY
srun --mpi=none --cpu-bind=cores --distribution=block:block bash -ceu '
    pair=$((SLURM_PROCID+1))
    read -r _ path _ score < <(sed -n "$((pair+1))p" "$out/manifest.tsv")
    awk "/Cpus_allowed_list:/{print \$2}" /proc/self/status > "$out/cpus-$pair.txt"
    "$work/trial" "$data/$path" "$pair" "$score" > "$out/pair-$pair.csv"
'
