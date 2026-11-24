#!/usr/bin/env python3
"""Run paired fresh-process trials of the same CUDA kernels with two allocation policies."""
import argparse
import csv
import hashlib
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path

def fixed_fasta(path, length, records):
    rng = random.Random(2026 + length)
    ref = rng.choices("ACGT", k=length)
    with path.open("w") as stream:
        for i in range(records):
            seq = list(ref)
            for j in rng.sample(range(length), round(0.02 * length)):
                seq[j] = rng.choice("ACGT".replace(seq[j], ""))
            stream.write(f">tip{i}\n{''.join(seq)}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("executable", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--cohort", type=Path, default=Path("data/phix174/cohort.fa"))
    args = p.parse_args()
    fields = ("job", "case", "input_sha256", "trial", "mode", "pairs", "seconds", "malloc_calls",
              "free_calls", "reuse_hits", "allocated_bytes", "peak_reserved_bytes", "peak_live_bytes",
              "malloc_seconds", "free_seconds", "distance_fnv64")
    with tempfile.TemporaryDirectory() as work, args.output.open("x", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        out.writeheader()
        cases = []
        for length, count in ((256, 16), (5000, 8), (20000, 8)):
            name = "fixed-" + str(length)
            path = Path(work) / (name + ".fa")
            fixed_fasta(path, length, count)
            cases.append((name, path))
        cases.append(("phix174", args.cohort))
        for case, path in cases:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            for trial in range(1, 4):
                modes = ("normal", "reuse") if trial % 2 else ("reuse", "normal")
                for mode in modes:
                    raw = subprocess.check_output([args.executable, path, mode], text=True)
                    result = json.loads(raw)
                    out.writerow(dict(result, job=os.environ.get("SLURM_JOB_ID", "local"),
                                      case=case, input_sha256=digest, trial=trial, mode=mode))
                    stream.flush()


if __name__ == "__main__":
    main()
