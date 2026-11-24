"""End-to-end NvPD OpenMP matrix comparison on a reproducible synthetic cohort."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("initial", type=Path, help="original intra-pair CLI")
    p.add_argument("candidate", type=Path, help="CLI with inter-pair scheduling")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--bind", action="store_true", help="pin each worker to one allocated core")
    args = p.parse_args()
    env = os.environ | {"OMP_NUM_THREADS": "8", "OMP_DYNAMIC": "FALSE",
                        "OMP_PROC_BIND": "FALSE", "OMP_MAX_ACTIVE_LEVELS": "4"}
    for name in ("OMP_WAIT_POLICY", "OMP_PLACES", "GOMP_CPU_AFFINITY"):
        env.pop(name, None)
    if args.bind:
        env.update(OMP_PROC_BIND="close", OMP_PLACES="cores")
    rng = random.Random(20261124)
    ancestor = "".join(rng.choices("ACGT", k=2048))
    sequences = ["".join(rng.choice("ACGT".replace(c,"")) if rng.random()<0.03 else c
                         for c in ancestor) + "A"*(i%7) for i in range(16)]
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp,"cohort.fa")
        path.write_text("".join(f">seq{i:02d}\n{s}\n" for i,s in enumerate(sequences)))
        input_sha = hashlib.sha256(path.read_bytes()).hexdigest()

        def run(binary, algorithm):
            start = time.perf_counter()
            out = subprocess.check_output([binary,"matrix","--algorithm",algorithm,"--threads","8",
                                           "--input",path], env=env, timeout=120)
            return time.perf_counter()-start, out

        _, oracle = run(args.candidate,"edlib")
        matrix_sha = hashlib.sha256(oracle).hexdigest()
        layouts = [("intra", args.initial), ("inter", args.candidate)]
        for trial in range(6):
            # Discard the first invocation of both binaries, including its startup.
            for order, (layout,binary) in enumerate(layouts if trial%2 else reversed(layouts),1):
                elapsed, matrix = run(binary,"nvpd-omp")
                if matrix != oracle:
                    raise ValueError(f"matrix differs from edlib for {layout}, trial {trial}")
                if trial:
                    rows.append(dict(trial=trial, order=order, layout=layout, threads=8, records=16,
                                     pairs=120, wall_seconds=elapsed, input_sha256=input_sha,
                                     matrix_sha256=matrix_sha))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    paired = {r["trial"]: {} for r in rows}
    for r in rows:
        paired[r["trial"]][r["layout"]] = r["wall_seconds"]
    print(f"Median paired speedup {statistics.median(v['intra']/v['inter'] for v in paired.values()):.3f}x")


if __name__ == "__main__":
    main()
