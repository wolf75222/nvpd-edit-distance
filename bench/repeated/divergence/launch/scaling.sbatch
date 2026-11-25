#!/bin/bash -l
#SBATCH --job-name=nvpd-grace-scaling
#SBATCH --partition=instant
#SBATCH --constraint=armgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
set -u
out=$(realpath -m "${1:?output directory}")
panel=$(realpath "${2:-data/human/chm13-panel}")
mkdir "$out"
work=$(mktemp -d "${TMPDIR:-/dev/shm}/chps-grace.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
bench="$PWD/bench/repeated/divergence/launch"
record="$PWD/results/human/chromosome-panel"

cmake -S . -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    -DCHPS_WARNINGS_AS_ERRORS=ON -DCMAKE_PROJECT_INCLUDE="$bench/scaling.cmake" >"$out/configure.log"
cmake --build "$work/build" --parallel 4 --target grace_scaling >"$out/build.log"

python3 - "$panel" "$record" "$out" <<'PY'
import csv
import hashlib
import sys
from pathlib import Path

panel, record, out = map(Path, sys.argv[1:])
with (record / "trials.csv").open() as stream:
    trials = list(csv.DictReader(stream))
with (record / "manifest.tsv").open() as stream:
    manifest = {row["case"]: row for row in csv.DictReader(stream, delimiter="\t")}
with (out / "inputs.tsv").open("w") as inputs, (out / "inputs.sha256").open("w") as hashes:
    for case in ("chr1", "chr2", "chr3", "chr4"):
        row = manifest[case]
        scores = {method: {int(r["distance"]) for r in trials if r["case"] == case and r["method"] == method}
                  for method in ("edlib", "wfa2")}
        if len(scores["wfa2"]) != 1 or scores["wfa2"] != scores["edlib"]:
            raise ValueError(f"Reference scores disagree for {case}")
        path = panel / (case + ".fa")
        digest = hashlib.sha256()
        with path.open("rb") as fasta:
            for chunk in iter(lambda: fasta.read(4 << 20), b""):
                digest.update(chunk)
        if digest.hexdigest() != row["input_sha256"]:
            raise ValueError(f"Different input for {case}")
        print(path, row["m"], row["n"], scores["wfa2"].pop(), sep="\t", file=inputs)
        print(digest.hexdigest(), path.name, file=hashes)
PY

{
    hostname -s
    printf 'job=%s\n' "$SLURM_JOB_ID"
    g++ --version | head -1
    lscpu | awk '/^Model name:|^Architecture:|^Thread\(s\) per core:|L2 cache:|L3 cache:/'
    sha256sum "$bench"/scaling.{cpp,cmake,sbatch} "$work/build/grace_scaling"
} >"$out/environment.txt"

srun --mpi=none --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
    python3 - "$work/build/grace_scaling" "$out" <<'PY'
import csv
import os
from pathlib import Path
import random
import subprocess
import sys

binary, out = Path(sys.argv[1]), Path(sys.argv[2])
physical = {}
for cpu in sorted(os.sched_getaffinity(0)):
    topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    key = tuple((topology / name).read_text().strip() for name in ("physical_package_id", "core_id"))
    physical.setdefault(key, cpu)
if len(physical) != 4:
    raise RuntimeError(f"Expected four physical cores in the Slurm step, found {physical}")
cpus = list(physical.values())
env = dict(os.environ, OMP_PLACES=",".join(f"{{{cpu}}}" for cpu in cpus),
           OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close", OMP_MAX_ACTIVE_LEVELS="1")
with (out / "environment.txt").open("a") as stream:
    print("physical_package/core:logical_cpu", physical, file=stream)
order = list(range(6))
random.Random(840841).shuffle(order)
with (out / "trials.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    for block, permutation in enumerate(order, 1):
        command = [str(binary), str(out / "inputs.tsv"), str(block), str(permutation)]
        with subprocess.Popen(command, stdout=subprocess.PIPE, text=True, env=env) as proc:
            rows = iter(csv.reader(proc.stdout))
            header = next(rows)
            if block == 1:
                writer.writerow(header)
            count = 0
            for row in rows:
                writer.writerow(row)
                stream.flush()
                count += 1
            if proc.wait() or count != 9:
                raise RuntimeError(f"Block {block} incomplete")
PY
