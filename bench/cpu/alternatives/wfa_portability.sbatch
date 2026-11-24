#!/bin/bash -l
#SBATCH --job-name=wfa-portability
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
tmp=$(mktemp -d "${TMPDIR:-/tmp}/wfa-ab.XXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
out="$PWD/results/alternatives/wfa-portability/$SLURM_JOB_ID"
mkdir -p "$out" "$tmp/source"
cp -R CMakeLists.txt src include cmake "$tmp/source/"
[[ ! -d deps ]] || ln -s "$PWD/deps" "$tmp/source/deps"
sha256sum CMakeLists.txt cmake/patches/wfa2-portability.cmake bench/quicked/warm.cpp \
    src/libraries/baselines.cpp src/libraries/wfa2/wfa_adapter.c > "$out/sources.sha256"
{ hostname -s; lscpu | grep -E 'Model name:|Thread.s. per core:';
  gcc --version | head -1; objdump --version | head -1; echo 'C/C++ flags -O3 -DNDEBUG -march=native'; } > "$out/environment.txt"
echo 'variant,object,part,sha256' > "$out/machine-code.csv"
for v in stock portable; do
    # Both builds use the identical source/build path, flags and library revision.
    mkdir "$tmp/$v"
    if [[ $v == stock ]]; then : > "$tmp/source/cmake/patches/wfa2-portability.cmake"
    else cp cmake/patches/wfa2-portability.cmake "$tmp/source/cmake/patches/"; fi
    cmake -S "$tmp/source" -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
        -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native > "$tmp/config-$v.log"
    cmake --build "$tmp/build" --parallel 4 --target editdist > "$tmp/build-$v.log"
    g++ -O3 -DNDEBUG -march=native -std=c++20 -fopenmp -DCHPS_BENCH_NATIVE -Iinclude \
        bench/quicked/warm.cpp src/io/fasta.cpp "$tmp/build/libdistances.a" "$tmp/build/libedlib_static.a" \
        "$tmp/build/libwfa_static.a" -lm -o "$tmp/$v/warm"
    for unit in wavefront_extend_kernels mm_allocator editdist warm; do
        case "$unit" in
            editdist) obj="$tmp/build/editdist" ;;
            warm) obj="$tmp/$v/warm" ;;
            *) obj=$(find "$tmp/build/CMakeFiles/wfa_static.dir" -name "$unit.c.o");;
        esac
        [[ -f $obj ]]
        sha=$(sha256sum "$obj" | cut -d' ' -f1)
        echo "$v,$unit,ELF,$sha" >> "$out/machine-code.csv"
        while read -r section; do
            objcopy -O binary --only-section="$section" "$obj" "$tmp/section"
            sha=$(sha256sum "$tmp/section" | cut -d' ' -f1)
            echo "$v,$unit,$section,$sha" >> "$out/machine-code.csv"
        done < <(objdump -h "$obj" | awk '$1 ~ /^[0-9]+$/ {name=$2} /CODE/ {print name}')
        if [[ $unit == wavefront_extend_kernels || $unit == mm_allocator ]]; then
            objdump -dr "$obj" | sed '/file format/d' > "$tmp/$v/$unit.dis"
            sha=$(sha256sum "$tmp/$v/$unit.dis" | cut -d' ' -f1)
            echo "$v,$unit,disassembly-and-relocations,$sha" >> "$out/machine-code.csv"
        fi
    done
    rm -r -- "$tmp/build"
done
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE RAYON_NUM_THREADS=1
export WFA_BENCH_OUT="$out"
srun --mpi=none --exact -N1 -n1 -c1 --cpu-bind=cores bash -c '
    taskset -pc $$ | sed "s/pid [0-9]*.s current affinity list: /bound core /" >> "$WFA_BENCH_OUT/environment.txt"
    exec python3 bench/cpu/alternatives/wfa_portability.py "$@"
' bash "$tmp" "${CHPS_WFA_INPUTS:-data/controls/quicked}" "$out"
