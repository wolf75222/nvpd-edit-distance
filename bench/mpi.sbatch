#!/bin/bash -l
#SBATCH --job-name=chps-mpi-scale
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
spack load /lni3xp2
pmix_prefix=$(spack location -i /zp25yhu)
export LD_LIBRARY_PATH="$pmix_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set -u
build=$(mktemp -d "${TMPDIR:-/tmp}/chps-mpi.XXXXXXXX")
trap 'rm -rf -- "$build"' EXIT
study=${1:-strong}
case "$study" in
    strong) out="results/mpi/$SLURM_JOB_ID" ;;
    weak) out="results/mpi/weak/$SLURM_JOB_ID" ;;
    *) echo 'expected strong or weak' >&2; exit 2 ;;
esac
mkdir -p "$out"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE
{
    date -u +%FT%TZ
    scontrol show hostnames "$SLURM_JOB_NODELIST"
    g++ --version | head -1
    mpirun --version | head -1
    cmake --version | head -1
    scontrol show job -o "$SLURM_JOB_ID" | tr ' ' '\n' | \
        grep -E '^(JobId|Partition|NumNodes|NumCPUs|NumTasks|CPUs/Task|NodeList|TimeLimit|ReqTRES|AllocTRES)='
} > "$out/environment.txt"
sha256sum src/mpi.cpp src/dp.cpp src/dp_internal.hpp src/fasta.cpp bench/mpi* > "$out/source.sha256"
cmake -S . -B "$build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native" \
    -DCHPS_MPI=ON -DCHPS_CUDA=OFF -DBUILD_TESTING=ON
cmake --build "$build" --target editdist-mpi editdist check_distances --parallel 12
ctest --test-dir "$build" -R '^distances$' --output-on-failure > "$out/native-tests.txt"
# Node-local builds need one copy on each node. Use a unique shared path
# just for the executable and remove it when the allocation ends.
shared=$(mktemp -d "$PWD/.mpi-bin.XXXXXXXX")
trap 'rm -rf -- "$build" "$shared"' EXIT
cp "$build/editdist-mpi" "$shared/"
export CHPS_SLURM_MPI=${CHPS_SLURM_MPI:-pmix}
export CHPS_MPI_DATA_TMP="$shared"
if [[ -n ${KOKKOS_SOURCE:-} ]]; then
    # A failed auxiliary build is reported separately. It must not suppress MPI timing.
    KOKKOS_THREADS=8 srun --exclusive --exact --mpi=none -N1 -n1 -c8 --cpu-bind=cores bash bench/kokkos-step.sh \
        cpu "results/kokkos/cpu-$SLURM_JOB_ID" > "$out/kokkos-step.log" 2>&1 || echo "Kokkos step failed, see its log"
fi
mapfile -t nodes < <(scontrol show hostnames "$SLURM_JOB_NODELIST")
export CHPS_MPI_HOSTS="${nodes[0]},${nodes[1]}"
printf '>a\nACG\n>b\nACTG\n' > "$shared/probe.fa"
probe=("$shared/editdist-mpi" pair --input "$shared/probe.fa")
check_probe() {
    "$@" "${probe[@]}" > "$out/probe.json" 2> "$out/probe.log" &&
        python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["distance"]==1; rows=[json.loads(x) for x in open(sys.argv[2]) if x.startswith("{")]; assert len(rows)==1 and rows[0]["ranks"]==2' "$out/probe.json" "$out/probe.log"
}
if check_probe timeout 45 srun --exclusive --exact --mpi="$CHPS_SLURM_MPI" -N2 -n2 --ntasks-per-node=1 -c1 --cpu-bind=cores; then
    export CHPS_MPI_LAUNCHER=srun
    smoke="srun --exclusive --exact --mpi=$CHPS_SLURM_MPI -N1 -n4 -c2 --cpu-bind=cores"
else
    mv "$out/probe.log" "$out/srun-probe.log"
    export CHPS_MPI_LAUNCHER=mpirun
    check_probe timeout 60 mpirun --host "${nodes[0]}:192,${nodes[1]}:192" -n 2 --map-by ppr:1:node --bind-to core
    smoke="mpirun --host ${nodes[0]}:192 -n 4 --map-by ppr:4:node:PE=2 --bind-to core"
fi
echo "launcher=$CHPS_MPI_LAUNCHER plugin=$CHPS_SLURM_MPI" >> "$out/environment.txt"
python3 bench/mpi_checks.py --binary "$shared/editdist-mpi" --oracle "$build/editdist" \
    --launcher "$smoke" > "$out/mpi-tests.txt"
if [[ $study == weak ]]; then
    python3 bench/mpi_weak.py --binary "$shared/editdist-mpi" --oracle "$build/editdist" --out "$out"
    exit
fi
cohort=()
if [[ -f data/phix174/cohort.fa ]]; then cohort=(--cohort data/phix174/cohort.fa); fi
python3 bench/mpi_benchmark.py --binary "$shared/editdist-mpi" --oracle "$build/editdist" \
    --out "$out" "${cohort[@]}"
