#!/bin/bash -l
#SBATCH --job-name=chps-1m-openmp
#SBATCH --nodes=2
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:45:00
#SBATCH --output=completion-1m-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
revision=f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-1m.XXXXXXXX")
shared=$(mktemp -d "$PWD/.chps-1m-bin.XXXXXXXX")
trap 'rm -rf -- "$work" "$shared"' EXIT
base=${BASELINE_DIR:-$work/initial}
if [[ -z ${BASELINE_DIR:-} ]]; then
    mkdir "$base"
    git archive "$revision" CMakeLists.txt include src tests | tar -xC "$base"
    if [[ -d deps ]]; then ln -s "$PWD/deps" "$base/deps"; fi
fi
export out="$PWD/results/romeo/completion-1m/jobs/$SLURM_JOB_ID"
mkdir -p "$out/raw"
root=$PWD
(cd "$base" && sha256sum -c "$root/results/romeo/sources.sha256") > "$out/source-check.txt"
cmake -S "$base" -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DCHPS_CUDA=OFF
cmake --build "$work/build" --parallel 8
ctest --test-dir "$work/build" --output-on-failure
export binary="$shared/editdist"
cp "$work/build/editdist" "$binary"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE
{
    echo "source_revision=$revision"
    echo "job=$SLURM_JOB_ID"
    scontrol show hostnames "$SLURM_JOB_NODELIST"
    lscpu | awk -F: '/^Model name:/{print}'
    g++ --version | head -1
    cmake --version | head -1
    awk '/^CXX_FLAGS =/{print}' "$work/build/CMakeFiles/distances.dir/flags.make"
    echo "binary_sha256=$(sha256sum "$binary" | cut -d' ' -f1)"
    echo "OMP_PLACES=$OMP_PLACES OMP_PROC_BIND=$OMP_PROC_BIND OMP_DYNAMIC=$OMP_DYNAMIC"
} > "$out/environment.txt"
for rep in 1 2 3 4 5; do
    export rep
    srun --mpi=none --nodes=2 --ntasks=6 --ntasks-per-node=3 --cpus-per-task=8 \
        --cpu-bind=cores --distribution=block:block bash -c '
        set -eu
        rates=(0.02 0.05 0.1)
        alg=diagonal
        if (( (rep + SLURM_NODEID) % 2 == 0 )); then alg=nvpd-omp; fi
        stem="$out/raw/rep-$rep-node-$SLURM_NODEID-task-$SLURM_LOCALID"
        awk "/^Cpus_allowed_list:/{print \$2}" /proc/self/status > "$stem.cpus"
        python3 bench/campaign/run.py --binary "$binary" --target cpu --profile full --lengths 1000000 \
            --rates "${rates[SLURM_LOCALID]}" --replicates "$rep" --algorithms "$alg" \
            --threads 8 --trials 1 --output "$stem.csv"
    '
    echo "replicate $rep complete $(date -u +%FT%TZ)"
done
srun --mpi=none --exact --nodes=1 --ntasks=3 --ntasks-per-node=3 --cpus-per-task=8 \
    --cpu-bind=cores --distribution=block:block bash -c '
    set -eu
    rates=(0.02 0.05 0.1)
    python3 bench/campaign/run.py --binary "$binary" --target cpu --profile full --lengths 1000000 \
        --rates "${rates[SLURM_LOCALID]}" --algorithms edlib wfa2 --trials 1 \
        --output "$out/raw/refs-$SLURM_LOCALID.csv"
'
python3 bench/cpu/openmp/million/collect_1m.py "$out" results/romeo
