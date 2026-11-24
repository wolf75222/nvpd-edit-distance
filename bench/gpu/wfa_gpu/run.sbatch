#!/bin/bash -l
#SBATCH --job-name=chps-wfagpu-affine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=24G
#SBATCH --constraint=armgpu
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
runner=bench/gpu/wfa_gpu
out=results/wfa-gpu/$SLURM_JOB_ID
mkdir -p "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-wfagpu.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
unpack() {
    local name=$1 path=$2 url=$3 hash=$4 archive="deps/$1.tar.gz"
    if [[ ! -f $archive ]]; then
        archive="$work/$name.tar.gz"
        curl -fsSL "$url" -o "$archive"
    fi
    echo "$hash  $archive" | sha256sum -c --status
    mkdir -p "$path"
    tar xf "$archive" -C "$path" --strip-components=1
}
gpu_rev=85acbcf4224811b3e52814b1bbe9364b5e2ba5ed
wfa_rev=49c255df126ee536fe92caff7a9f7c183ec3ff29
unpack wfa-gpu "$work/upstream" "https://codeload.github.com/Marco-Sola-Lab/WFA-GPU/tar.gz/$gpu_rev" e67e94177e36d4782c7c656f274522d445c507a4503d18f5d45484e8345018c6
unpack wfa-gpu-wfa2 "$work/upstream/external/WFA2-lib" "https://codeload.github.com/smarco/WFA2-lib/tar.gz/$wfa_rev" 14432a1e7e40f2c4df2354ab5e56cbf3761f97a6f5c5065cf62fe8538f1bc118
# Expose the upstream fallback counter in release builds; the threshold control exercises it.
sed -i 's/LOG_DEBUG("(Batch %d) %d\/%d alignemnts/LOG_INFO("(Batch %d) %d\/%d alignemnts/' "$work/upstream/lib/align.cu"
# The upstream CLI links with gcc. Resolve its indirect CUDA dependency from Spack.
cuda_root=$(dirname "$(dirname "$(readlink -f "$(command -v nvcc)")")")
cudart=$(find -L "$cuda_root" -type f -name 'libcudart.so' -print -quit)
[[ -n $cudart ]] || { echo "libcudart absent under $cuda_root" >&2; exit 1; }
cuda_lib=$(dirname "$cudart")
export LD_LIBRARY_PATH="$cuda_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
make -C "$work/upstream" "ARGS_ALIGNER=-Wall -Wno-unused-function -Lbuild/ -I. -Ilib/ -Wl,-rpath-link,$cuda_lib -Wl,-rpath,$cuda_lib" SM=90 COMPUTE=90 BUILD_TOOLS=0 BUILD_EXAMPLES=0
cmake -S . -B "$work/baselines" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 "-DCMAKE_PROJECT_INCLUDE=$PWD/$runner/targets.cmake"
cmake --build "$work/baselines" -j 4 --target wfagpu_baselines
{
    hostname -s
    nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv
    nvcc --version
    gcc --version | head -1
    echo "WFA-GPU $gpu_rev; WFA2 CPU fallback $wfa_rev"
    echo "edlib ec2310eda1841ab48c14cd3d866778a4f5eb1491; WFA2-lib bcf473a6561fa297934a80eaa7c04b4ee412360c"
} > "$out/environment.txt"
sha256sum "$runner"/{batch.cpp,targets.cmake,measure.py,run.sbatch} scripts/datasets/synthetic.py src/{cuda/cuda.cu,libraries/baselines.cpp,libraries/wfa2/wfa_adapter.c} > "$out/source.sha256"
if [[ ${1:-} == --sars ]]; then
    (($# == 2)) && ((SLURM_CPUS_PER_TASK >= 5)) || { echo 'require --sars FASTA and five cores' >&2; exit 2; }
    sha256sum "$runner/sars_cohort.py" results/ncbi/sars_content.tsv results/phylogeny/sars-edit/matrix.tsv >> "$out/source.sha256"
    srun --mpi=none --cpu-bind=cores python3 "$runner/sars_cohort.py" "$work" "$out" "$2"
elif (($#)); then
    sha256sum "$runner/chromosome.py" results/mpi/human/{windows,scores}.tsv >> "$out/source.sha256"
    srun --mpi=none --cpu-bind=cores python3 "$runner/chromosome.py" "$work" "$out" "$1"
else
    srun --mpi=none --cpu-bind=cores python3 "$runner/measure.py" "$work" "$out"
fi
