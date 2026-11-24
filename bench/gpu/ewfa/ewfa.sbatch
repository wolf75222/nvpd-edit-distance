#!/bin/bash -l
#SBATCH --job-name=chps-ewfa-2022
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
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
out="results/ewfa/$SLURM_JOB_ID"
mkdir -p "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/chps-ewfa.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
rev=8b7817cc1623b297aace606aeb0b72dca64d39a8
archive="${EWFA_ARCHIVE:-$work/ewfa2022.tar.gz}"
if [[ ! -f $archive ]]; then
    curl -fsSL "https://codeload.github.com/quim0/eWFA-GPU/tar.gz/$rev" -o "$archive"
fi
echo "0d2510090270d707b1ad25dd84ae8aa4269040102b3c4718d2e19312111c4f38  $archive" | sha256sum --check --status
edlib=deps/edlib.tar.gz
if [[ ! -f $edlib ]]; then
    edlib="$work/edlib.tar.gz"
    curl -fsSL https://codeload.github.com/Martinsos/edlib/tar.gz/ec2310eda1841ab48c14cd3d866778a4f5eb1491 -o "$edlib"
fi
echo "ca2fc2fe74ad7d8342094b15e201ead2d63c8f5014a4376d2ef623ee8a23d7ff  $edlib" | sha256sum --check --status
mkdir "$work/upstream" "$work/edlib"
tar xf "$archive" -C "$work/upstream" --strip-components=1
tar xf "$edlib" -C "$work/edlib" --strip-components=1
make -C "$work/upstream" SM=90 COMPUTE=90 all
cp "$work/upstream/wfe.aligner" "$work/original"
patch -d "$work/upstream" -p1 < bench/gpu/ewfa/ewfa_validation.patch
make -C "$work/upstream" SM=90 COMPUTE=90 all
cp "$work/upstream/wfe.aligner" "$work/raw"
nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -Iinclude -Isrc/cpu \
    -I"$work/edlib/edlib/include" bench/gpu/ewfa/ewfa_batch.cpp src/cuda/cuda.cu src/cpu/dp.cpp \
    "$work/edlib/edlib/src/edlib.cpp" -o "$work/baselines"
{
    hostname -s
    nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv
    nvcc --version
    echo "eWFA-GPU $rev"
    echo "edlib ec2310eda1841ab48c14cd3d866778a4f5eb1491"
} > "$out/environment.txt"
sha256sum src/cpu/dp.cpp src/cpu/dp_internal.hpp src/cuda/cuda.cu bench/gpu/ewfa/ewfa.sbatch bench/gpu/ewfa/ewfa_batch.cpp bench/gpu/ewfa/ewfa_validation.patch bench/gpu/ewfa/ewfa_measure.py scripts/datasets/synthetic.py > "$out/source.sha256"
python3 scripts/datasets/synthetic.py --lengths 5000 --rates 0.005 0.02 0.05 --seed 2026 --pairs 5 --out "$work/synthetic"
python3 bench/gpu/ewfa/ewfa_measure.py --work "$work" --out "$out"
