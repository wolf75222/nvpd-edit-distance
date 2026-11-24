#!/bin/bash -l
#SBATCH --job-name=chr21-exact-path
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
tmp=$(mktemp -d "/project/$SLURM_JOB_ACCOUNT/$USER/quicked-human.XXXXXXXX")
mkdir -p "$tmp/tmp"
export TMPDIR="$tmp/tmp"
trap 'rm -rf -- "$tmp"' EXIT
out=results/alternatives/human/$SLURM_JOB_ID
mkdir -p "$out" "$tmp/astar" "$tmp/cargo"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=false RAYON_NUM_THREADS=1
# Same pinned compiler, Cargo.lock and source revision as bench/quicked/run.sbatch.
revision=$(git get-tar-commit-id < <(gzip -dc deps/astar-src.tar.gz))
[[ $revision == bf2e14e0cbc3a9a03600dcda0641d7f89e401e63 ]]
tar -xzf deps/astar-src.tar.gz -C "$tmp/astar"
export RUSTUP_HOME="$tmp/rustup" CARGO_HOME="$tmp/cargo" PATH="$tmp/cargo/bin:$PATH"
tar --warning=no-unknown-keyword -xzf deps/cargo-cache.tar.gz -C "$CARGO_HOME"
(cd "$tmp/astar" && echo 'ab82655822401810e4a29e77b9e02eda9c67284f7706e6e1d224de9b2cd7e12c  Cargo.lock' | sha256sum -c -)
curl -fsSL --retry 2 https://static.rust-lang.org/dist/2026-05-21/channel-rust-nightly.toml -o "$tmp/channel.toml"
echo "6c7ed921781dbc26497f3a06f05bd02307aa8a70952c85dd5f043334051ec5ab  $tmp/channel.toml" | sha256sum -c -
curl -fsSL --retry 2 --connect-timeout 12 --max-time 120 https://sh.rustup.rs -o "$tmp/rustup-init.sh"
sh "$tmp/rustup-init.sh" -y --no-modify-path --profile minimal --default-toolchain nightly-2026-05-21
(cd "$tmp/astar" && cargo build --release --locked -j 4 -p astarpa-c)
build="$tmp/build"
cmake -S . -B "$build" -DCHPS_QUICKED=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_C_FLAGS=-march=native -DCMAKE_CXX_FLAGS=-march=native
cmake --build "$build" --parallel 4 --target distances quicked_static
printf '{global: chps_quicked_certified; local: *;};\n' >"$tmp/symbols.map"
g++ -shared -Wl,-Bsymbolic -Wl,--version-script="$tmp/symbols.map" -Wl,--whole-archive "$build/libquicked_static.a" \
    -Wl,--no-whole-archive "$build/libedlib_static.a" -lm -o "$tmp/libcertified.so"
g++ -std=c++20 -O3 -DNDEBUG -march=native -fopenmp -Iinclude -I"$tmp/astar/astarpa-c" \
    bench/quicked/human/trial.cpp src/io/fasta.cpp "$build/libdistances.a" "$build/libedlib_static.a" "$build/libwfa_static.a" \
    -L"$tmp" -lcertified -L"$tmp/astar/target/release" -lastarpa_c -Wl,-rpath,"$tmp:$tmp/astar/target/release" -lm -o "$tmp/trial"
{ hostname -s; gcc --version | head -1; rustc --version; lscpu | grep -E 'Model name|Thread.s. per core'; } >"$out/environment.txt"
(cd "$tmp" && sha256sum trial libcertified.so astar/target/release/libastarpa_c.so) >"$out/binaries.sha256"
sha256sum bench/quicked/human/{run.sbatch,trial.cpp,measure.py} src/libraries/quicked/quicked_adapter.c src/libraries/wfa2/wfa_adapter.c src/libraries/baselines.cpp src/io/fasta.cpp >"$out/sources.sha256"
srun --mpi=none --exact -n1 -c4 --cpu-bind=cores python3 bench/quicked/human/measure.py "$tmp/trial" data/human "$out"
