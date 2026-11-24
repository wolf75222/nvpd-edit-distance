#!/bin/bash -l
#SBATCH --job-name=nvpd-pages
#SBATCH --constraint=x64cpu
#SBATCH --partition=instant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%x-%j.log
set -eo pipefail
romeo_load_x64cpu_env >/dev/null
spack load /q22mihh
set -u
panel=${1:?provide the directory containing the CHM13 panel}
perf=${PERF:-/project/r260105/rmdraux/tools/perf/bin/perf}
out="results/profile/huge-pages/$SLURM_JOB_ID"
mkdir -p "$out"
tmp=$(mktemp -d "${TMPDIR:-$PWD}/pages.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
cmake -S . -B "$tmp/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native > "$tmp/config.log"
cmake --build "$tmp/build" --parallel "$SLURM_CPUS_PER_TASK" --target distances > "$tmp/build.log"
g++ -std=c++20 -O3 -march=native -fopenmp -Iinclude bench/cpu/huge-pages/trial.cpp src/io/fasta.cpp "$tmp/build"/{libdistances,libedlib_static,libwfa_static}.a -lm -o "$tmp/trial"
export OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=true OMP_DYNAMIC=false
{
    hostname -s
    uname -sr
    g++ --version | head -1
    "$perf" --version
    lscpu | awk '/^Model name:|^L[123].*cache:/'
    for f in /sys/kernel/mm/transparent_hugepage/{enabled,defrag,hpage_pmd_size} /proc/sys/kernel/perf_event_paranoid; do echo "$f $(cat "$f")"; done
    getconf LEVEL1_DCACHE_LINESIZE
} > "$out/environment.txt"
"$perf" list --details > "$tmp/events" 2>/dev/null
grep -E -A3 '^  ls_(l1_d_tlb_miss.all_l2_miss|(any_fills_from_sys|dmnd_fills_from_sys|hw_pf_dc_fills).dram_io_(near|far))' "$tmp/events" > "$out/events.txt"
{ "$perf" stat -a -e power/energy-pkg/k -- sleep .2 || true
  "$perf" stat -a -M umc_mem_read_bandwidth -- sleep .2 || true
  for f in /sys/class/{powercap/*/energy_uj,hwmon/hwmon*/energy*_input}; do [[ ! -r "$f" ]] || echo "$f"; done
} > "$out/socket-counters.txt" 2>&1
(cd "$panel" && sha256sum chr1.fa chr18.fa chr20.fa) > "$out/inputs.sha256"
sha256sum src/cpu/split.cpp bench/cpu/huge-pages/{trial.cpp,run.sbatch} > "$out/sources.sha256"
echo 'case,block,solver,pages,m,n,distance,mapped_bytes,huge_bytes,setup_s,calls,core_active,cycles,instructions,l1d_miss,dtlb_walks,core_s,mem_active,all_fills,dem_fills,prefetch_fills,mem_s' > "$out/trials.csv"
for pair in chr1:176 chr18:5 chr20:30; do
    name=${pair%:*}
    srun -N1 -n1 -c1 --cpu-bind=cores "$tmp/trial" "$name" "$panel/$name.fa" "${pair#*:}" >> "$out/trials.csv"
done
echo 'block,pages,queries,distance,setup_s,total_s,mapped_bytes,huge_bytes' > "$out/reuse.csv"
srun -N1 -n1 -c1 --cpu-bind=cores "$tmp/trial" chr1 "$panel/chr1.fa" 176 reuse >> "$out/reuse.csv"
