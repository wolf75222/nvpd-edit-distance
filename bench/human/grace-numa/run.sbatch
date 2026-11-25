#!/bin/bash -l
#SBATCH --job-name=nvpd-interchip
#SBATCH --partition=instant
#SBATCH --constraint=armgpu
#SBATCH --nodelist=romeo-a002
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-socket=1
#SBATCH --cores-per-socket=1
#SBATCH --cpus-per-task=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
set -u
root=${CHPS_ROOT:-$PWD}
out=$(realpath -m "${1:?output directory}")
manifest=$(realpath "${2:?manifest containing paths to the four chromosome FASTA files}")
study="$root/bench/human/grace-numa"
mkdir "$out"
work=$(mktemp -d "${TMPDIR:-/tmp}/grace-numa.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
{
    hostname -s
    date -Iseconds
    scontrol show job -dd "$SLURM_JOB_ID"
    lscpu -e=CPU,NODE,SOCKET,CORE | awk 'NR==1 || $1==0 || $1==72'
    cat /proc/sys/kernel/numa_balancing
    g++ --version | head -1
    sha256sum "$study"/{numa.cpp,numa.cmake,measure.py,run.sbatch}
} >"$out/environment.txt"
# A FASTA may have moved since the first campaign. Check each byte against its original digest.
python3 - "$manifest" "$root/results/human/grace-scaling/inputs.sha256" "$out" <<'PY'
import hashlib
import sys
from pathlib import Path
manifest, checksums, out = map(Path, sys.argv[1:])
expected = {name: digest for digest, name in (line.split() for line in checksums.read_text().splitlines())}
rows = [line.split("\t") for line in manifest.read_text().splitlines()]
if [Path(row[0]).name for row in rows] != [f"chr{i}.fa" for i in range(1, 5)]:
    raise ValueError("Expected original CHM13 chromosomes 1 to 4 in order")
for name, *_ in rows:
    path = Path(name)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 << 20), b""):
            digest.update(chunk)
    if digest.hexdigest() != expected[path.name]:
        raise ValueError(f"FASTA digest mismatch: {path.name}")
(out / "inputs.sha256").write_bytes(checksums.read_bytes())
PY
cmake -S "$root" -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCHPS_WARNINGS_AS_ERRORS=ON -DCMAKE_PROJECT_INCLUDE="$study/numa.cmake" >"$out/configure.log"
cmake --build "$work/build" --parallel 2 --target grace_numa >"$out/build.log"
sha256sum "$work/build/grace_numa" >>"$out/environment.txt"
taskset -c 0,72 python3 "$study/measure.py" "$work/build/grace_numa" "$manifest" "$out"
