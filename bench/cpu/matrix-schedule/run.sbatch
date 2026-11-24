#!/bin/bash -l
#SBATCH --job-name=chps-rows
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=matrix-schedule-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
input=${1:-data/sars/cohort.fa}
matrix=results/phylogeny/sars-edit/matrix.tsv
run_dir="results/phylogeny/sars-schedule/$SLURM_JOB_ID"
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-schedule.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
mkdir -p "$work/source" "$run_dir"
cp -R CMakeLists.txt cmake include src "$work/source/"
[[ ! -d deps ]] || ln -s "$PWD/deps" "$work/source/deps"
source_cli="$work/source/src/cli/main.cpp"
if [[ -n ${MATRIX_BASELINE_MAIN:-} ]]; then
    cp "$MATRIX_BASELINE_MAIN" "$source_cli"
else
    git show d9f394c:src/cli/main.cpp > "$source_cli"
fi
cmake -S "$work/source" -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native
cmake --build "$work/build" --target editdist --parallel 8
cp "$work/build/editdist" "$work/rows"
rows_source_sha=$(sha256sum "$source_cli" | cut -d' ' -f1)
python3 - "$source_cli" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
old = 'const std::size_t groups = pair_parallel ?'
new = ('const bool library_groups = across && n > 1 && '
       '(o.algorithm == "edlib" || o.algorithm == "wfa2" || o.algorithm == "biwfa");\n    '
       'const std::size_t groups = (pair_parallel || library_groups) ?')
source = p.read_text()
if source.count(old) != 1:
    raise ValueError("unexpected initial matrix scheduler")
p.write_text(source.replace(old, new))
PY
cmake --build "$work/build" --target editdist --parallel 8
cp "$work/build/editdist" "$work/groups"
groups_source_sha=$(sha256sum "$source_cli" | cut -d' ' -f1)
{
    hostname -s
    g++ --version | head -1
    lscpu | awk -F: '/^Model name:/{print}'
    echo 'Release -O3 -DNDEBUG -march=native'
    sha256sum "$input" "$matrix"
    (cd "$work" && sha256sum rows groups)
    echo "$rows_source_sha  rows:src/cli/main.cpp"
    echo "$groups_source_sha  groups:src/cli/main.cpp"
} > "$run_dir/environment.txt"
export OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
srun --mpi=none --exact -N1 -n1 -c32 --hint=nomultithread --cpu-bind=cores \
     python3 bench/cpu/matrix-schedule/compare.py "$work/rows" "$work/groups" "$input" "$matrix" "$run_dir/trials.csv"
