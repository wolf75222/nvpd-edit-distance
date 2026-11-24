#!/bin/bash -l
#SBATCH --job-name=chps-cpu-profile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --constraint=x64cpu
#SBATCH --output=profile-%j.log
set -euo pipefail
phase=${1:-baseline}
out=${PROFILE_OUTPUT:-results/profile/cpu/$SLURM_JOB_ID}
mkdir -p "$out"
out=$(realpath "$out")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/chps-profile.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
base=${BASELINE_DIR:-$tmp/initial}
if [[ -z ${BASELINE_DIR:-} ]]; then
    mkdir "$base"
    git archive f08b4f6 -- include/chps/dp.hpp src/dp.cpp src/dp_internal.hpp | tar -xC "$base"
fi
compiler=${CXX:-g++}
verify_source() {
    local file=$1 expected=$2
    if [[ $(sha256sum "$file" | cut -d' ' -f1) != "$expected" ]]; then
        echo "Source digest mismatch: $(basename "$file")" >&2
        exit 1
    fi
}
verify_source "$base/src/dp.cpp" 95a75e9e3dc1d1540aafa73ba05313269fdfbcd3740453d169a3c041446d63c6
verify_source "$base/src/dp_internal.hpp" 06711e4b010cc9254e8b0f6f8e14fc7ff93afc522388e5b3110954b9e9e0eddb
verify_source "$base/include/chps/dp.hpp" 5e423542698c662a7f7e5fffe3409860e8b516e8dc7b9ec8fa58cdfe52872b9d
flags=(-std=c++20 -O3 -g -march=native -fopenmp -fno-omit-frame-pointer)
build_variant() {
    local name=$1 src=$2
    "$compiler" "${flags[@]}" -I"$base/include" -I"$base/src" -fopt-info-vec-all="$tmp/$name.vec" \
        "$src" bench/profile_cpu_driver.cpp -o "$tmp/$name"
    awk '/dp.cpp:|v[0-9]+.cpp:/{if ($0 ~ /loop vectorized|couldn.t vectorize loop|control flow in loop/) {sub(/^.*\//,"");print}}' \
        "$tmp/$name.vec" | head -c 100000 > "$out/vector-$name.txt"
    sha256sum "$src" | cut -d' ' -f1 > "$out/source-$name.sha256"
}
build_variant baseline "$base/src/dp.cpp"
if [[ $phase != baseline ]]; then
    v1=${V1_SRC:-$tmp/v1.cpp}
    if [[ -z ${V1_SRC:-} ]]; then git show 4724f68:src/dp.cpp > "$v1"; fi
    verify_source "$v1" eb8dd38c31d714c0477bc5479a8301b532c6bbfc71cf924738ca9ec57dc70dc1
    build_variant v1 "$v1"
    if [[ $phase == v2 ]]; then
        v2=${V2_SRC:-$tmp/v2.cpp}
        if [[ -z ${V2_SRC:-} ]]; then git show experiment/compact-index:src/dp.cpp > "$v2"; fi
        verify_source "$v2" 57a9cfe89822c07eca43af863864cb70ab13c329d3ac0b238d937df10119470a
        build_variant v2 "$v2"
    fi
fi
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
{
    date -u +%FT%TZ
    hostname -s
    echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
    "$compiler" --version | head -n1
    echo "flags=${flags[*]}; OMP_PLACES=$OMP_PLACES OMP_PROC_BIND=$OMP_PROC_BIND OMP_WAIT_POLICY=$OMP_WAIT_POLICY"
    lscpu | awk -F: '/Model name:|L1d cache:|L2 cache:|L3 cache:|NUMA node\(s\):/{print}'
    for f in /sys/devices/system/cpu/cpu0/cache/index*/{level,type,size,shared_cpu_list}; do
        echo "${f#/sys/devices/system/cpu/cpu0/cache/} $(cat "$f")"
    done
    echo "perf_event_paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)"
} > "$out/environment.txt"
lscpu --extended=CPU,CORE,SOCKET,NODE,CACHE > "$out/topology.txt"
srun --ntasks=1 --cpu-bind=cores --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    python3 bench/profile_cpu.py "$phase" "$tmp" "$out"
