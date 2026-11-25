#!/usr/bin/env python3
"""Time the same NvPD calls in balanced process order, outside the profiler."""
import argparse
import csv
import hashlib
import io
import math
from pathlib import Path
import random
import subprocess
from statistics import mean


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("binaries", type=Path, nargs="+")
    parser.add_argument("--blocks", type=int, default=10)
    parser.add_argument("--lengths", nargs="+", type=int, default=(16384, 32768))
    parser.add_argument("--calls", type=int, default=1, help="timed calls after one warm-up per process")
    parser.add_argument("--threads", nargs="+", type=int, default=(1, 4))
    args = parser.parse_args()
    if args.blocks < 1 or args.calls < 1:
        raise ValueError("Positive block and call counts required")
    with (args.inputs / "inputs.csv").open(newline="") as stream:
        cases = [case for case in csv.DictReader(stream) if int(case["length"]) in args.lengths]
    if {int(case["length"]) for case in cases} != set(args.lengths):
        raise ValueError("Missing requested length in manifest")
    if len({p.name for p in args.binaries}) != len(args.binaries):
        raise ValueError("Binary names must identify the variants")
    for case in cases:
        path = args.inputs / f"{case['cohort']}_{case['length']}.fa"
        if hashlib.sha256(path.read_bytes()).hexdigest() != case["fasta_sha256"]:
            raise ValueError(f"Input differs from oracle: {path}")
    order = args.binaries.copy()
    rng = random.Random(256832)
    rng.shuffle(order)
    conditions = [(case, threads) for case in cases for threads in args.threads]
    with args.out.open("x", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        shape = ("columns",) if "columns" in cases[0] else ()
        writer.writerow(("cohort", "length", *shape, "threads", "block", "position", "variant", "score",
                         "warmup_seconds", "seconds"))
        for block in range(args.blocks):
            rng.shuffle(conditions)
            for case, threads in conditions:
                # Each condition has the same rotation across blocks regardless of its order in the block.
                for pos in range(len(order)):
                    binary = order[(pos + block) % len(order)]
                    path = args.inputs / f"{case['cohort']}_{case['length']}.fa"
                    columns = [case["columns"]] if shape else []
                    text = subprocess.check_output(
                        [binary.resolve(), path, case["length"], str(threads), str(args.calls + 1), *columns], text=True)
                    records = list(csv.DictReader(io.StringIO(text)))
                    if len(records) != args.calls + 1 or any(r["score"] != case["score"] or
                            not (math.isfinite(float(r["seconds"])) and float(r["seconds"]) > 0)
                            for r in records):
                        raise ValueError(f"Invalid result for {binary}: {records}")
                    writer.writerow((case["cohort"], case["length"], *columns, threads, block+1, pos+1,
                                     binary.name, case["score"], records[0]["seconds"], mean(float(r["seconds"]) for r in records[1:])))
                    stream.flush()
            print(f"Completed block {block+1}", flush=True)


if __name__ == "__main__":
    main()
