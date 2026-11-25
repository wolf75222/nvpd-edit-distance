#!/bin/bash -l
#SBATCH --job-name=chps-cqa-final
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=x64cpu
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=maqao-%j.log
set -euo pipefail
base=$(realpath "${1:?reference source tree}")
prod=$(realpath "${2:?production source tree}")
inputs=$(realpath "${3:?prepared pair FASTAs}")
compare=$(realpath "${4:?A/B results directory}")
variants=$(realpath "${5:?prepared experimental source directory}")
mkdir "${6:?new output directory}"
out=$(realpath "$6")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-cqa-final.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=4
maqao=${MAQAO:-maqao}
common=(-O3 -g -fopenmp -std=c++20 -Wall -Wextra -Werror -Wpedantic)
for mode in native generic; do
    flags=("${common[@]}")
    [[ $mode == generic ]] || flags+=(-march=native)
    mkdir "$work/$mode"
    for variant in baseline production; do
        root=$base
        [[ $variant == baseline ]] || root=$prod
        g++ "${flags[@]}" -I"$root/include" -I"$root/src/cpu" "$root/src/cpu/dp.cpp" \
            "$root/src/io/fasta.cpp" bench/cpu/maqao/kernel.cpp -o "$work/$mode/$variant"
    done
    g++ "${flags[@]}" -I"$prod/include" -I"$prod/src/cpu" "$prod/src/cpu/dp.cpp" \
        "$variants/check_variant.cpp" -o "$work/check-$mode"
    srun --mpi=none --cpu-bind=cores "$work/check-$mode" > "$out/validation-$mode.txt"
done
{
    hostname -s
    g++ --version | head -1
    "$maqao" --version | head -1
    echo "COMMON_FLAGS=${common[*]}; NATIVE_FLAGS=-march=native"
    echo "OMP_PROC_BIND=$OMP_PROC_BIND OMP_PLACES=$OMP_PLACES OMP_WAIT_POLICY=$OMP_WAIT_POLICY"
    srun --mpi=none --cpu-bind=cores bash -c 'grep Cpus_allowed_list /proc/self/status'
    sha256sum "$base/src/cpu/dp.cpp" "$prod/src/cpu/"{dp.cpp,nvpd_avx512.hpp,dp_internal.hpp}
} > "$out/environment.txt"
# Same FA, repeat counts and LProf configuration as the baseline profile in compare.sbatch.
while IFS=, read -r threads selected iterations; do
    [[ $threads != threads ]] || continue
    xp="$work/production-$threads"
    export OMP_NUM_THREADS=$threads
    srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html \
        --lprof-params='sampling-rate=high' -xp="$xp" -- "$work/native/production" \
        "$inputs/near_32768.fa" 32768 "$threads" "$iterations" < /dev/null > "$out/oneview-production-$threads.log" 2>&1
    tar -C "$xp" -czf "$out/production-$threads.tar.gz" .
done < "$compare/selection.txt"
export OMP_NUM_THREADS=4
for mode in native generic; do
    bins=("$work/$mode/baseline" "$work/$mode/production")
    srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$out/medium-$mode.csv" \
        "${bins[@]}" --lengths 32768 --blocks 8
done
bins=("$work/native/baseline" "$work/native/production")
srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$out/short-native.csv" \
    "${bins[@]}" --lengths 32 128 512 1024 4096 --blocks 8 --calls 12
srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$out/long-native.csv" \
    "${bins[@]}" --lengths 100000 --blocks 8
