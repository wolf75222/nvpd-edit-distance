#!/bin/bash -l
#SBATCH --job-name=chps-chromosomes
#SBATCH --partition=short
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=panel-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
data=$(realpath "${1:-data/human/chm13-panel}")
out="${PANEL_OUTPUT:-results/human/chromosome-panel/$SLURM_JOB_ID}"
mkdir -p "$out"
out=$(realpath "$out")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-chromosomes.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
if [[ -n ${PANEL_INITIAL:-} ]]; then
    cp "$PANEL_INITIAL" "$work/initial.cpp"
else
    git show a728923:src/io/fasta.cpp > "$work/initial.cpp"
fi
echo "28242b98ce45a329bdc94a54e7fec7cdde53a030ab407a171408466710e528f0  $work/initial.cpp" | sha256sum -c --status
sed '/#pragma omp simd reduction/d' src/io/fasta.cpp > "$work/scalar.cpp"
bench=bench/human/chromosome-panel
cmake -S . -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-march=native \
    -DCMAKE_C_FLAGS=-march=native -DCMAKE_PROJECT_INCLUDE="$PWD/$bench/targets.cmake" \
    -DPANEL_INITIAL="$work/initial.cpp" -DPANEL_SCALAR="$work/scalar.cpp"
cmake --build "$work/build" --parallel 8 --target panel_initial panel_scalar panel_simd
export OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
cp "$work/build/parser.vec" "$out/vectorization.txt"
{
    hostname -s
    g++ --version | head -1
    lscpu | awk '/^Model name:|L1d cache:|L2 cache:|L3 cache:|NUMA node\(s\):/'
    if [[ -r /proc/sys/kernel/perf_event_paranoid ]]; then
        cat /proc/sys/kernel/perf_event_paranoid || true
    fi
    sha256sum "$bench"/{run.sbatch,measure.py,trial.cpp,targets.cmake,counters.hpp} src/io/fasta.cpp
    sha256sum "$work/build"/panel_* | awk '{n=split($2,a,"/"); print $1,a[n]}'
} > "$out/environment.txt"
cp "$data/manifest.tsv" "$out/manifest.tsv"
srun --cpu-bind=cores --distribution=block:block python3 "$bench/measure.py" "$data" "$work/build" "$out"
