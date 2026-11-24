#!/bin/bash -l
#SBATCH --job-name=nvpd-nsys
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --constraint=armgpu
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
spack load /phhwuh3
set -u
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
tmp=$(mktemp -d "${TMPDIR:-/tmp}/nvpd-trace.XXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
out="$PWD/results/repeated/gpu-trace/$SLURM_JOB_ID"
mkdir -p "$out"
pair="$PWD/data/synthetic/pairs/n5000_r0.02_p1_f7364a035b8b7d3d.fa"
echo "4861b0646f36eb6f71939a9425cb757a6e85a9f792b0c35c85d16705e9c59201  $pair" | sha256sum -c --status
source="$tmp/source"
if [[ -n ${CHPS_TRACE_SOURCE:-} ]]; then source="$CHPS_TRACE_SOURCE"
else bash bench/repeated/trace/gpu_trace_sources.sh "$source"; fi
(cd "$source" && sha256sum -c sources.sha256)
cp "$source/sources.sha256" "$out/"
{ nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv; nvcc --version; nsys --version;
  echo "job=$SLURM_JOB_ID cpus=$SLURM_CPUS_PER_TASK gpu=1:sm_90"; } > "$out/environment.txt"
for v in baseline candidate; do
    nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -I"$source/include" -I"$source/src" \
        "$source/bench/profile_gpu.cu" "$source/$v.cu" "$source/src/fasta.cpp" "$source/src/dp.cpp" -o "$tmp/$v"
done
mkdir -p "$out/raw"
for block in {1..7}; do
    order=(baseline candidate)
    if (( block % 2 == 0 )); then order=(candidate baseline); fi
    for v in "${order[@]}"; do
        name="$out/raw/$block-$v"
        report="$tmp/$block-$v"
        nsys profile --trace=cuda --sample=none --cpuctxsw=none --capture-range=cudaProfilerApi \
            --capture-range-end=stop --output="$report" \
            bash -c 'output=$1; shift; exec "$@" > "$output"' bench "$name.log" \
            "$tmp/$v" "$v" "$pair" 1 "$block" 0 1 99 > "$tmp/nsys.log"
        nsys stats --report cuda_gpu_kern_sum,cuda_api_sum,cuda_gpu_mem_time_sum --format csv \
            --output "$name" "$report.nsys-rep" > /dev/null
    done
done
python3 bench/repeated/trace/gpu_trace.py "$out/raw" "$out"
python3 -m zipfile -c "$out/nsys.zip" "$out/raw"
rm -r -- "$out/raw"
