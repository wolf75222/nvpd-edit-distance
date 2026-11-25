#!/bin/bash -l
#SBATCH --job-name=nvpd-energy
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=12G
#SBATCH --time=00:35:00
#SBATCH --output=energy-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
set -u
inputs=$(realpath "${1:?provide bench/kokkos/human/inputs.py FASTA directory}")
base=bench/gpu/energy
out=results/energy/gh200/$SLURM_JOB_ID
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-energy.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
export OMP_NUM_THREADS=1
flags=(-O3 -std=c++20 -arch=sm_90 -DCHPS_CUDA_PACKED_MATCHES=0 -Iinclude -Isrc/cuda -Isrc/cpu)
nvcc "${flags[@]}" -c "$base/trial.cu" -o "$tmp/trial.o"
g++ -O3 -std=c++20 -fopenmp -Iinclude -c src/cpu/dp.cpp -o "$tmp/dp.o"
g++ -O3 -std=c++20 -Iinclude -c src/io/fasta.cpp -o "$tmp/fasta.o"
for rows in 4 8 16 32; do
    nvcc "${flags[@]}" -DCHPS_CUDA_TILE_ROWS="$rows" src/cuda/cuda.cu "$tmp/"{trial,dp,fasta}.o \
        -Xcompiler=-fopenmp,-pthread -l:libnvidia-ml.so.1 -o "$tmp/trial-$rows"
done
nvcc "${flags[@]}" bench/gpu/roofline/calibrate.cu -o "$tmp/calibrate"
{
    hostname -s
    nvidia-smi --query-gpu=name,uuid,pci.bus_id,driver_version,power.limit,memory.total --format=csv
    nvidia-smi --query-compute-apps=gpu_uuid,process_name --format=csv
    nvcc --version |tail -2
    ncu --version |tail -1
    nsys --version
    grep -i 'ProfilingAdminOnly' /proc/driver/nvidia/params || true
    sha256sum "$base"/{trial.cu,run.sbatch,measure.py} src/cuda/{cuda.cu,local_tiles.cuh,packed_match.cuh} bench/gpu/roofline/calibrate.cu
} > "$out/environment.txt"
echo 'method,block,steps,updates,requested_bytes,seconds' > "$out/calibration.csv"
for block in 1 2 3 4 5; do "$tmp/calibrate" "$block" >> "$out/calibration.csv"; done
srun --mpi=none --cpu-bind=cores python3 "$base/measure.py" "$tmp" "$inputs" \
    results/kokkos/human/references.csv "$out/trials.csv"
case_path="$inputs/collinear-65536.fa"
set +e
ncu --kernel-name 'regex:.*row_tile.*' --launch-count 1 --metrics \
    dram__bytes_read.sum,dram__bytes_write.sum,smsp__inst_executed.sum,gpu__time_duration.sum \
    --csv --page raw --log-file "$out/ncu.txt" "$tmp/trial-16" "$case_path" 51 0
ncu_rc=$?
set -e
if (( ncu_rc != 0 )) && ! grep -q 'ERR_NVGPUCTRPERM' "$out/ncu.txt"; then exit "$ncu_rc"; fi
for rows in 4 8 16 32; do
    nsys profile --trace=cuda --sample=none --cpuctxsw=none -o "$tmp/activity-$rows" \
        "$tmp/trial-$rows" "$case_path" 51 0 >"$tmp/nsys-$rows.log" 2>&1
    nsys export -t sqlite -o "$tmp/activity-$rows.sqlite" "$tmp/activity-$rows.nsys-rep" > /dev/null
done
python3 "$base/activity.py" "$tmp" "$out/activity.csv"
