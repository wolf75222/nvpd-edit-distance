#!/usr/bin/env python3
"""Sample resident host memory of the focused alignments in a separate run."""
import argparse
import csv
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from focus import GROUPS, METHODS, prepare, sha256
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "memory"))
from grace_memory import one


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("original", "candidate", "editdist", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--groups", choices=GROUPS, nargs="+", default=("100k-2pct", "1m-2pct"))
    parser.add_argument("--blocks", type=int, default=3)
    args = parser.parse_args()
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < 2 or args.blocks < 1 or len(set(args.groups)) != len(args.groups):
        parser.error("Use at least two cores, a positive block count and distinct groups")
    binaries = {m: (args.original if m == "original" else args.candidate).resolve() for m in METHODS}
    hashes = {m: sha256(binary) for m, binary in binaries.items()}
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = "1"
    with TemporaryDirectory(prefix="nvpd-rss-") as temp:
        directory = Path(temp)
        cases = prepare(directory, args.groups, args.editdist.resolve())
        os.sched_setaffinity(0, {cpus[1]})
        with (args.output / "smaps.csv").open("x", newline="") as stream:
            writer = None
            for block in range(1, args.blocks + 1):
                for i, case in enumerate(cases):
                    shape = {"case": case["path"], "m": case["length"], "n": case["mutant_length"],
                             "distance": str(case["distance"])}
                    for pos, offset in enumerate((0, 1, 3, 2), 1):
                        method = METHODS[(i + block - 1 + offset) % 4]
                        algorithm = method if method in ("edlib", "wfa2") else "nvpd-cuda"
                        measured = one(binaries[method], directory / case["path"], shape, algorithm, cpus[0])
                        row = {"job": os.environ.get("SLURM_JOB_ID", "local"), "group": case["group"],
                               "rep": case["rep"], "sha256": case["sha256"], "distance": case["distance"],
                               "block": block, "order": pos, "method": method, "binary_sha256": hashes[method],
                               "worker_cpu": cpus[0], "observer_cpu": cpus[1], **measured}
                        if writer is None:
                            writer = csv.DictWriter(stream, fieldnames=row, lineterminator="\n")
                            writer.writeheader()
                        writer.writerow(row)
                    stream.flush()
                print(f"Finished RSS block {block}/{args.blocks}", flush=True)


if __name__ == "__main__":
    main()
