#!/bin/bash -l
#SBATCH --job-name=chps-cqa-ab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=x64cpu
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=maqao-%j.log
set -euo pipefail
root=$(realpath "${1:?original source tree}")
variants=$(realpath "${2:?directory with centered.cpp, blocked.cpp and blocked_rows.hpp}")
inputs=$(realpath "${3:?prepared pair FASTAs}")
mkdir "${4:?new output directory}"
out=$(realpath "$4")
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-cqa-ab.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
flags=(-O3 -g -march=native -fopenmp -std=c++20 -Wall -Wextra -Werror -I"$root/include" -I"$root/src/cpu")
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE OMP_NUM_THREADS=4
maqao=${MAQAO:-"$HOME/.local/bin/maqao"}
for v in baseline wide centered centered-wide blocked; do
    cpp="$root/src/cpu/dp.cpp"; extra=()
    case "$v" in
        wide) extra=(-mprefer-vector-width=512);;
        centered) cpp="$variants/centered.cpp";;
        centered-wide) cpp="$variants/centered.cpp"; extra=(-mprefer-vector-width=512);;
        blocked) cpp="$variants/blocked.cpp"; extra=(-Wno-maybe-uninitialized -Wno-uninitialized);;
    esac
    g++ "${flags[@]}" "${extra[@]}" "$cpp" "$root/src/io/fasta.cpp" bench/cpu/maqao/kernel.cpp -o "$work/$v"
done
g++ "${flags[@]}" -Wno-uninitialized -Wno-maybe-uninitialized "$variants/blocked.cpp" "$variants/check_variant.cpp" -o "$work/check"
srun --mpi=none --cpu-bind=cores "$work/check" > "$out/validation.txt"
{
    hostname -s
    g++ --version | head -1
    "$maqao" --version | head -1
    echo "flags=${flags[*]}"
    echo "OMP_PROC_BIND=$OMP_PROC_BIND OMP_PLACES=$OMP_PLACES OMP_WAIT_POLICY=$OMP_WAIT_POLICY"
    srun --cpu-bind=cores bash -c 'grep Cpus_allowed_list /proc/self/status'
    sha256sum "$root/src/cpu/"{dp.cpp,dp_internal.hpp} "$variants/"{centered.cpp,blocked.cpp,blocked_rows.hpp}
} > "$out/environment.txt"
bins=("$work"/{baseline,wide,centered,centered-wide,blocked})
srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$out/trials.csv" "${bins[@]}"
# Select on the 32k inputs, then validate on the longer, unprofiled extent.
python3 bench/cpu/maqao/select_variant.py "$out/trials.csv" > "$out/selection.txt"
mapfile -t candidates < <(tail -n +2 "$out/selection.txt" | cut -d, -f2 | sort -u)
held=("$work/baseline")
for v in "${candidates[@]}"; do [[ $v == baseline ]] || held+=("$work/$v"); done
blocks=8
(( ${#held[@]} != 3 )) || blocks=9
srun --mpi=none --cpu-bind=cores python3 bench/cpu/maqao/compare.py "$inputs" "$out/long.csv" "${held[@]}" --lengths 100000 --blocks "$blocks"
# Compare the original and best hot loops on the same CPU, each profiled for at least 15 s (target 22 s).
while IFS=, read -r threads selected iterations; do
    [[ $threads != threads ]] || continue
    export OMP_NUM_THREADS=$threads
    for v in baseline "$selected"; do
        key="$v-$threads"
        [[ ! -e "$out/$key.tar.gz" ]] || continue
        xp="$work/profile-$key"
        srun --mpi=none --cpu-bind=cores "$maqao" oneview -R1 --output-format=html \
            --lprof-params='sampling-rate=high' -xp="$xp" -- "$work/$v" "$inputs/near_32768.fa" \
            32768 "$threads" "$iterations" < /dev/null > "$out/oneview-$key.log" 2>&1
        tar -C "$xp" -czf "$out/$key.tar.gz" .
    done
done < "$out/selection.txt"
