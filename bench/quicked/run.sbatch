#!/bin/bash -l
#SBATCH --job-name=quicked-ab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --constraint=x64cpu
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE RAYON_NUM_THREADS=1
tmp=$(mktemp -d "/project/$SLURM_JOB_ACCOUNT/$USER/quicked-ab.XXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
out=results/alternatives/romeo/$SLURM_JOB_ID
mkdir -p "$out"
test -s inputs/quicked/cases.csv || { echo "Prepare the FASTA files with bench/quicked/cases.py inputs/quicked" >&2; exit 2; }
test -f deps/astar-src.tar.gz || { echo "Stage the pinned A*PA archive (see results/alternatives/README.md)." >&2; exit 2; }
cmake -S . -B "$tmp/cxx" -DCMAKE_BUILD_TYPE=Release -DCHPS_QUICKED=ON \
    -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native
cmake --build "$tmp/cxx" -j8 --target editdist-quicked editdist check_quicked_widen measure-quicked measure-quicked-stock measure-quicked-native
ctest --test-dir "$tmp/cxx" --output-on-failure -R '^quicked'

build_astar() (
    set -e
    mkdir "$tmp/astar" "$tmp/cargo"
    revision=$(git get-tar-commit-id < <(gzip -dc deps/astar-src.tar.gz))
    [[ $revision == bf2e14e0cbc3a9a03600dcda0641d7f89e401e63 ]]
    tar -xzf deps/astar-src.tar.gz -C "$tmp/astar"
    export RUSTUP_HOME="$tmp/rustup" CARGO_HOME="$tmp/cargo" PATH="$tmp/cargo/bin:$PATH"
    if [[ -f deps/cargo-cache.tar.gz ]]; then tar -xzf deps/cargo-cache.tar.gz -C "$CARGO_HOME"; fi
    (cd "$tmp/astar" && echo "ab82655822401810e4a29e77b9e02eda9c67284f7706e6e1d224de9b2cd7e12c  Cargo.lock" | sha256sum -c -)
    curl -fsSL --retry 2 https://static.rust-lang.org/dist/2026-05-21/channel-rust-nightly.toml -o "$tmp/channel.toml"
    echo "6c7ed921781dbc26497f3a06f05bd02307aa8a70952c85dd5f043334051ec5ab  $tmp/channel.toml" | sha256sum -c -
    # The pinned crate uses nightly features. Install only inside the job.
    curl -fsSL --retry 2 --connect-timeout 12 --max-time 120 https://sh.rustup.rs -o "$tmp/rustup-init.sh"
    sh "$tmp/rustup-init.sh" -y --no-modify-path --profile minimal --default-toolchain nightly-2026-05-21
    (cd "$tmp/astar" && cargo build --release --locked -j 6 -p astarpa-c)
    g++ -std=c++20 -O3 -march=native -DCHPS_BENCH_ASTAR -Iinclude -I"$tmp/astar/astarpa-c" \
        bench/quicked/warm.cpp src/fasta.cpp -L"$tmp/astar/target/release" -lastarpa_c \
        -Wl,-rpath,"$tmp/astar/target/release" -o "$tmp/astar-trial"
    rustc --version > "$out/rust-version.txt"
)
build_astar
# Keep originals for stock, including its architecture-specific failures.
cp -R inputs/quicked "$tmp/inputs"
g++ -std=c++20 -O3 -march=native -Iinclude -I"$tmp/cxx/_deps/edlib-src/edlib/include" \
    bench/quicked/find_witness.cpp "$tmp/cxx/CMakeFiles/measure-quicked-stock.dir/_deps/quicked-src/quicked/src/bpm_banded.c.o" \
    "$tmp/cxx/libquicked_static.a" "$tmp/cxx/libedlib_static.a" -lm -o "$tmp/find-witness"
"$tmp/find-witness" "$tmp/inputs" > "$tmp/found.csv"
while IFS=, read -r name length exact score; do
    mkdir -p "$out/witnesses"
    cp "$tmp/inputs/$name.fa" "$out/witnesses/"
    sha=$(sha256sum "$tmp/inputs/$name.fa" | cut -d' ' -f1)
    echo "$name,$length,$length,$exact,$sha" >> "$tmp/inputs/cases.csv"
done < "$tmp/found.csv"
python3 bench/quicked/diagnose.py "$tmp/cxx" "$tmp/inputs" "$out/trace.csv"
export BENCH_OUT="$out" BENCH_BIN="$tmp/cxx" BENCH_INPUT="$tmp/inputs" BENCH_ALLOC="$SLURM_CPUS_PER_TASK"
srun --mpi=none --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores bash -c '
    { echo "job $SLURM_JOB_ID"; echo "allocated_cpus $BENCH_ALLOC"; echo "active_threads 1";
      hostname -s; lscpu | grep -E "Model name:|Thread.s. per core:";
      taskset -pc $$; g++ --version | head -1; echo "flags -O3 -DNDEBUG -march=native"; } > "$BENCH_OUT/environment.txt"
    python3 bench/quicked/measure.py --inputs "$BENCH_INPUT" --bin "$BENCH_BIN" --out "$BENCH_OUT" "$@"
' bash --astar "$tmp/astar-trial"
