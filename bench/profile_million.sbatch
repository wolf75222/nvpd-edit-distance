#!/bin/bash -l
#SBATCH --job-name=chps-million-ab
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
set -u
export work=$(mktemp -d "${TMPDIR:-/tmp}/chps-million.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export out="results/profile/million/$SLURM_JOB_ID"
mkdir -p "$out"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
refs=(f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83 4724f68bf7dbfb83b133422541b165b5c7124da1)
variants=(initial simd)
digests=(95a75e9e3dc1d1540aafa73ba05313269fdfbcd3740453d169a3c041446d63c6 eb8dd38c31d714c0477bc5479a8301b532c6bbfc71cf924738ca9ec57dc70dc1)
for i in 0 1; do
    v=${variants[i]}
    if [[ -n ${CHPS_MILLION_VERSIONS:-} ]]; then
        cp "$CHPS_MILLION_VERSIONS/$v.cpp" "$work/$v.cpp"
    else
        git show "${refs[i]}:src/dp.cpp" > "$work/$v.cpp"
    fi
    echo "${digests[i]}  $work/$v.cpp" | sha256sum --check --status
    g++ -std=c++20 -O3 -DNDEBUG -march=native -fopenmp -Iinclude -Isrc bench/million_pair.cpp \
        "$work/$v.cpp" src/fasta.cpp -o "$work/$v"
done
{
    hostname -s
    g++ --version | head -1
    echo 'flags=-std=c++20 -O3 -DNDEBUG -march=native -fopenmp'
    echo "initial=${refs[0]} simd=${refs[1]}"
    sha256sum bench/million_pair.cpp bench/profile_million.sbatch src/{dp_internal.hpp,fasta.cpp}
} > "$out/environment.txt"
srun --mpi=none --exact -N1 -n3 -c8 --hint=nomultithread --cpu-bind=cores --distribution=block:block bash -c '
    set -eu
    rates=(0.02 0.05 0.1)
    rate=${rates[SLURM_LOCALID]}
    paths=(data/synthetic/pairs/n1000000_r${rate}_p1_*.fa)
    [[ ${#paths[@]} == 1 ]]
    path=${paths[0]}
    sha=$(sha256sum "$path" | cut -d" " -f1)
    variants=(initial simd)
    if (( SLURM_LOCALID%2 )); then variants=(simd initial); fi
    for pos in 0 1; do
        v=${variants[pos]}
        echo "start $rate $v $(date -u +%FT%TZ)" >&2
        values=$("$work/$v" "$path")
        echo "$SLURM_JOB_ID,$rate,$v,$((pos+1)),$sha,$values" >> "$work/group-$SLURM_LOCALID.csv"
    done'
echo 'job,rate,variant,order,input_sha256,m,n,distance,seconds,cpu_ids' > "$out/trials.csv"
cat "$work"/group-*.csv >> "$out/trials.csv"
python3 - "$out/trials.csv" <<'PY'
import csv, glob, sys
from pathlib import Path
rows=list(csv.DictReader(open(sys.argv[1])))
if len(rows)!=6: raise ValueError("missing timings")
oracle={}
for name in glob.glob("results/romeo/cpu-*.csv"):
    for x in csv.DictReader(open(name)):
        if x["length"]=="1000000" and x["rep"]=="1" and x["algorithm"] in ("edlib","wfa2") and x["status"]=="ok":
            oracle[float(x["rate"]),x["algorithm"]]=(x["input_sha256"],x["m"],x["n"],x["distance"])
occupied=set()
for rate in (.02,.05,.1):
    pair=[x for x in rows if float(x["rate"])==rate]
    if {x["variant"] for x in pair}!={"initial","simd"} or len({x["cpu_ids"] for x in pair})!=1:
        raise ValueError("paired layout mismatch")
    for x in pair:
        got=tuple(x[k] for k in ("input_sha256","m","n","distance"))
        if any(oracle[rate,v]!=got for v in ("edlib","wfa2")): raise ValueError("oracle discrepancy")
    physical=set()
    for cpu in pair[0]["cpu_ids"].split(";"):
        top=Path("/sys/devices/system/cpu/cpu"+cpu+"/topology")
        physical.add(tuple((top/key).read_text().strip() for key in ("physical_package_id","core_id")))
    if len(physical)!=8 or not physical.isdisjoint(occupied): raise ValueError("overlapping physical cores")
    occupied.update(physical)
print("6 exact results, scores and input digests agree with edlib and WFA2, 3 disjoint eight-core groups.")
PY
