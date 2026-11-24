#!/bin/bash -l
#SBATCH --job-name=chps-gpu-profile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1
#SBATCH --mem=8G
#SBATCH --constraint=armgpu
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /dvle736
spack load /phhwuh3
set -u
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
build=$(mktemp -d "${TMPDIR:-/tmp}/profile-gpu.XXXXXXXX")
trap 'rm -rf -- "$build"' EXIT
out="results/profile/gpu/$SLURM_JOB_ID"
mkdir -p "$out"
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv > "$out/device.txt"
nvcc --version >> "$out/device.txt"
nsys --version >> "$out/device.txt"
mode=${1:-profile}
pairs=data/synthetic/pairs
input=("$pairs/n5000_r0.02_p1_f7364a035b8b7d3d.fa" "$pairs/n20000_r0.02_p1_ed8d69e5ad60d1be.fa"
       "$pairs/n100000_r0.02_p1_c794ea53bcdb62d2.fa" "$pairs/n1000000_r0.02_p1_98b33c3d8b831451.fa")
expected=(99 396 1983 19801)
if [[ $mode == holdout ]]; then
    input=("$pairs/n5000_r0.1_p1_86a9646a353e3afe.fa" "$pairs/n20000_r0.1_p1_62ed7614335ca2f3.fa"
           "$pairs/n100000_r0.1_p1_79f31facc7f5111c.fa" "$pairs/n1000000_r0.1_p1_ee6350608d80f33c.fa")
    expected=(483 1908 9589 95720)
fi
export OMP_DISPLAY_ENV=FALSE
variants=(baseline)
ref=f08b4f6fcdf8
if [[ -z ${GPU_BASELINE:-} ]]; then
    git show "$ref:src/cuda.cu" > "$build/original.cu"
    GPU_BASELINE="$build/original.cu"
fi
sources=("$GPU_BASELINE")
echo "0c85e68539b93dd82deaaee744d3b72d2e76e44f906e78a1319d0b00dabf0f13  $GPU_BASELINE" \
    | sha256sum --check --status
if [[ $mode == compare || $mode == counters || $mode == holdout ]]; then
    variants+=(candidate)
    source_root=$PWD
    if [[ -n ${GPU_CANDIDATE_REV:-} ]]; then
        source_root="$build/source"
        mkdir "$source_root"
        cp -R CMakeLists.txt src include tests "$source_root/"
        ln -s "$PWD/bench" "$source_root/bench"
        if [[ -d deps ]]; then ln -s "$PWD/deps" "$source_root/deps"; fi
        git show "$GPU_CANDIDATE_REV:src/cuda.cu" > "$source_root/src/cuda.cu"
    fi
    sources+=("$source_root/src/cuda.cu")
    cmake -S "$source_root" -B "$build/cmake" -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON \
        -DCMAKE_CUDA_ARCHITECTURES=90 -DBUILD_TESTING=ON > "$build/config.log"
    cmake --build "$build/cmake" --parallel "$SLURM_CPUS_PER_TASK" > "$build/build.log"
    ctest --test-dir "$build/cmake" --output-on-failure | sed "s|$build/cmake|build|g" > "$out/validation.txt"
    "$build/cmake/editdist" pair --algorithm nvpd-cuda --input "${input[0]}" >> "$out/validation.txt"
    if [[ $mode == compare && ${GPU_INTERMEDIATE:-1} == 1 ]]; then
        python3 scripts/synthetic.py --out "$build/intermediate" --lengths 250000 500000 --rates 0.02 --pairs 1 --seed 2026
        for f in "$build"/intermediate/pairs/*.fa; do
            input+=("$f")
            expected+=("$("$build/cmake/editdist" pair --algorithm edlib --input "$f" \
                | python3 -c 'import json,sys; print(json.load(sys.stdin)["distance"])')")
        done
    elif [[ $mode == holdout ]]; then
        python3 bench/cohort.py --binary "$build/cmake/editdist" --fasta data/phix174/cohort.fa \
            --threads "$SLURM_CPUS_PER_TASK" --verify-with nvpd-cuda --output "$build/cohort"
        cp "$build/cohort/timings.csv" "$out/cohort.csv"
    fi
elif [[ $mode == profile ]]; then
    input=("${input[0]}")
else
    echo "expected profile, compare, counters or holdout" >&2
    exit 2
fi
sha256sum "${input[@]}" | sed -E 's|  .*/|  |' > "$out/inputs.sha256"
for i in "${!variants[@]}"; do
    # Both versions use the same host sources, compiler flags and input FASTAs.
    sha256sum "${sources[i]}" | awk -v name="${variants[i]}" '{print $1 "  " name}' >> "$out/variants.sha256"
    nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp -Iinclude -Isrc \
        bench/profile_gpu.cu "${sources[i]}" src/fasta.cpp src/dp.cpp -o "$build/${variants[i]}"
done
if [[ $mode == counters || ${GPU_COUNTERS:-0} == 1 ]]; then
    ncu --version >> "$out/device.txt"
    for v in "${variants[@]}"; do
        kernel=next_row
        if [[ $v == candidate ]]; then kernel=row_tile; fi
        # Profile one representative kernel from a warm 5k call on the visible GPU.
        if ! ncu --profile-from-start off --kernel-name "regex:$kernel" \
            --launch-count 1 --metrics 'gpu__time_duration.sum,sm__throughput.avg.pct_of_peak_sustained_elapsed,dram__bytes_read.sum,dram__bytes_write.sum' \
            --csv --page raw --log-file "$out/$v-counters.csv" \
            "$build/$v" "$v" "${input[0]}" 1 0 0 1 "${expected[0]}"; then
            if ! grep -q ERR_NVGPUCTRPERM "$out/$v-counters.csv"; then exit 1; fi
        fi
        if grep -q ERR_NVGPUCTRPERM "$out/$v-counters.csv"; then
            grep '^==ERROR==' "$out/$v-counters.csv" > "$out/$v-counters.txt"
            rm -- "$out/$v-counters.csv"
        fi
    done
    if [[ $mode == counters ]]; then exit 0; fi
fi
echo "variant,input,m,n,mode,run,trial,trace,seconds,distance,device_bytes" > "$out/trials.csv"
for i in "${!input[@]}"; do
    for run in 1 2 3; do
        order=("${variants[@]}")
        if ((run % 2 == 0 && ${#order[@]} == 2)); then order=(candidate baseline); fi
        for v in "${order[@]}"; do
            "$build/$v" "$v" "${input[i]}" 3 "$run" 0 0 "${expected[i]}" >> "$out/trials.csv"
        done
    done
done
if [[ $mode == holdout ]]; then exit 0; fi
for v in "${variants[@]}"; do
    # Captured runs are excluded from A/B wall timing.
    nsys profile --trace=cuda --sample=none --cpuctxsw=none --capture-range=cudaProfilerApi \
        --capture-range-end=stop --output="$build/$v" \
        "$build/$v" "$v" "${input[0]}" 1 0 0 1 "${expected[0]}" > "$build/$v.log"
    nsys stats --report cuda_gpu_kern_sum,cuda_api_sum,cuda_gpu_mem_time_sum --format csv \
        --output "$out/$v" "$build/$v.nsys-rep" >/dev/null
done
