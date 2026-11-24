#!/bin/bash -l
#SBATCH --job-name=chps-matrix-omp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
export work=$(mktemp -d "${TMPDIR:-/tmp}/chps-matrix.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export out="$PWD/results/profile/matrix-omp/$SLURM_JOB_ID"
mkdir -p "$out" "$work/source"
cp -R CMakeLists.txt include src "$work/source/"
[[ ! -d deps ]] || ln -s "$PWD/deps" "$work/source/deps"
cmake -S "$work/source" -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native
cmake --build "$work/build" --target editdist --parallel 8
cp "$work/build/editdist" "$work/inter"
rev=0cc2454b8c4778ddd8599e02d4f6e237c5e02650
if [[ -n ${MATRIX_INITIAL_MAIN:-} ]]; then
    cp "$MATRIX_INITIAL_MAIN" "$work/source/src/main.cpp"
else
    git show "$rev:src/main.cpp" > "$work/source/src/main.cpp"
fi
echo "0cfca4ac15f54df03cf1c7f966ac5a2fc788b83b6167547879f658dfa2b3a98e  $work/source/src/main.cpp" | sha256sum --check --status
cmake --build "$work/build" --target editdist --parallel 8
cp "$work/build/editdist" "$work/intra"
{
    hostname -s
    g++ --version | head -1
    lscpu | awk -F: '/^Model name:/{print}'
    echo 'Release -O3 -DNDEBUG -march=native'
    echo "original CLI revision $rev"
    sha256sum src/{main,dp,fasta}.cpp bench/{matrix_omp.py,matrix_omp.sbatch}
    (cd "$work" && sha256sum intra inter)
} > "$out/environment.txt"
srun --mpi=none --exact -N1 -n1 -c8 --hint=nomultithread --cpu-bind=cores bash -c '
    set -eu
    export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=4
    printf ">a\nACGTA\n>b\nAGTA\n" > "$work/probe.fa"
    for v in intra inter; do
        OMP_DISPLAY_AFFINITY=TRUE OMP_AFFINITY_FORMAT="%n %A" "$work/$v" pair --algorithm nvpd-omp \
            --threads 8 --input "$work/probe.fa" > /dev/null 2> "$work/$v-placement"
        sort -nu "$work/$v-placement" > "$out/$v-placement.txt"
    done
    cmp "$out/intra-placement.txt" "$out/inter-placement.txt"
    python3 bench/matrix_omp.py "$work/intra" "$work/inter" --bind --out "$out/trials.csv"'
python3 - "$out" <<'PY'
import csv,sys
from pathlib import Path
out=Path(sys.argv[1]); hosts=set(); workers=set()
for line in (out/"inter-placement.txt").read_text().splitlines():
    worker,cpu=map(int,line.split())
    top=Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    hosts.add(tuple((top/key).read_text().strip() for key in ("physical_package_id","core_id")))
    workers.add(worker)
if workers!=set(range(8)) or len(hosts)!=8: raise ValueError("unexpected placement")
rows=list(csv.DictReader((out/"trials.csv").open()))
if len(rows)!=10 or any(r["input_sha256"]!="e9070acc92776d65d26c9e5a80cdb4eec567997e5ca52ed8296e992c1ad05df1"
                       or r["matrix_sha256"]!="b1e413eb9070d56bca289459bce97bd9a5f32c55ec5ac99bb1ba43a0f9379ca8" for r in rows):
    raise ValueError("different synthetic input or edlib matrix")
print("Ten matrices agree with edlib on the original synthetic cohort; eight identical, distinct physical cores.")
PY
