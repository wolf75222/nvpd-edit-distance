#!/bin/bash -l
#SBATCH --job-name=chps-route-maqao
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=route-maqao-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
root=$(pwd -P)
fasta=$(realpath "${1:-build/hybrid-matrix/viral/surveillance.fa}")
reference=$(realpath "${2:-results/phylogeny/sars-edit/matrix.tsv}")
out=${3:-"build/maqao/route-$SLURM_JOB_ID"}
mkdir -p "$(dirname "$out")"
mkdir "$out"
out=$(realpath "$out")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-route-maqao.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
cc=${CC:-/usr/bin/gcc}
cxx=${CXX:-/usr/bin/g++}
flags="-O3 -g -DNDEBUG -march=native"
cmake -S "$root" -B "$work/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$cc" -DCMAKE_CXX_COMPILER="$cxx" -DCMAKE_C_FLAGS_RELEASE="$flags" \
    -DCMAKE_CXX_FLAGS_RELEASE="$flags" -DCHPS_WARNINGS_AS_ERRORS=ON -DBUILD_TESTING=OFF
cmake --build "$work/build" --parallel "$SLURM_CPUS_PER_TASK" --target editdist-route
binary="$work/build/editdist-route"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=FALSE
{
    printf 'job=%s\n' "$SLURM_JOB_ID"
    hostname -s
    lscpu | awk '/^(CPU\(s\)|Model name|Architecture):/'
    "$cxx" --version | head -1
    cmake --version | head -1
    printf 'flags=%s\n' "$flags"
    printf 'perf_event_paranoid='
    cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || printf 'unavailable\n'
    srun --mpi=none --cpu-bind=cores bash -c 'grep Cpus_allowed_list /proc/self/status'
    sha256sum "$binary" "$fasta" "$reference" src/hybrid/{matrix,main}.cpp src/hybrid/sketch.hpp
} > "$out/environment.txt"
profile() {
    local maqao=$1 destination=$2
    local xp="$work/profile" matrix="$work/profiled.tsv"
    "$maqao" --version > "$destination/version.txt"
    head -n 1 "$destination/version.txt" >> "$out/environment.txt"
    srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html,txt,csv \
        --lprof-params='sampling-rate=high' -xp="$xp" -- "$binary" --input "$fasta" \
        --policy hybrid-counts --threads "$SLURM_CPUS_PER_TASK" --warmups 1 --repeats 6 \
        --matrix "$matrix" --decisions "$destination/decisions.tsv" > "$destination/oneview.log" 2>&1
    cmp "$reference" "$matrix"
    (cd "$work" && sha256sum profiled.tsv) > "$destination/matrix.sha256"
    tar -C "$xp" -czf "$destination/profile.tar.gz" .
}

if (( $# < 4 )); then
    profile "${MAQAO:-maqao}" "$out"
else
    first_decisions=
    for version in "${@:4}"; do
        destination="$out/$version"
        mkdir "$destination"
        prefix=$(spack location -i "maqao@$version")
        profile "$prefix/bin/maqao" "$destination"
        if [[ -n $first_decisions ]]; then
            cmp "$first_decisions" "$destination/decisions.tsv"
        else
            first_decisions="$destination/decisions.tsv"
        fi
        mv "$work/profile" "$work/profile-$version"
    done
fi
