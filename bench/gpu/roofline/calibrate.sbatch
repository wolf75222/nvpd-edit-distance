#!/bin/bash -l
#SBATCH --job-name=nvpd-roof-cal
#SBATCH --constraint=armgpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=roof-cal-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
set -u
out=${ROOF_OUT:-results/roofline/gpu/$SLURM_JOB_ID}
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/roof-cal.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
src=bench/gpu/roofline/calibrate.cu
nvcc -O3 -std=c++20 -arch=sm_90 "$src" -o "$tmp/calibrate"
{
    hostname -s
    nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total,clocks.max.sm,clocks.max.memory --format=csv
    nvcc --version | tail -2
    ncu --version | tail -1
    nsys --version
    grep -i 'ProfilingAdminOnly' /proc/driver/nvidia/params || true
    sha256sum "$src" bench/gpu/roofline/calibrate.sbatch
} > "$out/environment.txt"
printf '%s\n' 'method,block,steps,updates,requested_bytes,seconds' > "$out/calibration.csv"
for block in 1 2 3 4 5; do "$tmp/calibrate" "$block" >>"$out/calibration.csv"; done
# Standard unprivileged collection. A driver-denied result is recorded verbatim.
set +e
ncu --kernel-name 'regex:stream_updates' --launch-count 1 --metrics \
    dram__bytes_read.sum,dram__bytes_write.sum,smsp__inst_executed.sum,gpu__time_duration.sum \
    --csv --page raw --log-file "$out/ncu-calibration.txt" "$tmp/calibrate" 0 >"$tmp/ncu.out" 2>&1
rc=$?
set -e
if ! grep -q 'ERR_NVGPUCTRPERM' "$out/ncu-calibration.txt" && (( rc != 0 )); then
    cat "$tmp/ncu.out"; exit "$rc"
fi
cuobjdump --dump-sass "$tmp/calibrate" | awk '/Function : .*register_updates/ {show=1;next} /Function : / {show=0} show' |
    sed -nE 's/.*\/\*[0-9a-fA-F]+\*\/ *(@!?P[0-9]+ +)?([A-Za-z0-9_.]+).*/\2/p' |
    sort | uniq -c >"$out/register-sass-counts.txt"
