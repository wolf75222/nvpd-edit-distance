#!/usr/bin/env python3
"""Compare production CUDA matrix commands from the baseline and reusable-workspace builds."""
import argparse
import csv
import hashlib
import json
import os
import random
import subprocess
import tempfile
import time
from pathlib import Path


def matrix(binary, fasta, backend):
    cmd = [binary, "matrix", "--algorithm", backend, "--input", fasta, "--threads", "4"]
    start = time.monotonic()
    result = subprocess.check_output(cmd)
    return time.monotonic() - start, result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("workspace", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--cohort", type=Path, default=Path("data/phix174/cohort.fa"))
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fasta = Path(temp) / "inputs.fa"
        rng = random.Random(2026)
        seqs = ["".join(rng.choices("ACGTRN", k=n)) for n in (32, 2048, 48, 4096)]
        fasta.write_text("".join(f">s{i}\n{x}\n" for i, x in enumerate(seqs)))
        ref = matrix(args.workspace, fasta, "edlib")[1]
        if not all(matrix(b, fasta, "nvpd-cuda")[1] == ref for b in (args.baseline, args.workspace)):
            raise ValueError("changing-width matrix mismatch")
        for a, b in [("", "N"), ("N" * 64 + "R", "A" * 516 + "R"), (seqs[0], seqs[2])]:
            fasta.write_text(f">a\n{a}\n>b\n{b}\n")
            cmd = ["pair", "--input", fasta, "--trace", "--algorithm"]
            oracle = json.loads(subprocess.check_output([args.workspace, *cmd, "seq"]))
            for exe in (args.baseline, args.workspace):
                val = json.loads(subprocess.check_output([exe, *cmd, "nvpd-cuda"]))
                if val["distance"] != oracle["distance"] or val["cigar"] != oracle["cigar"]:
                    raise ValueError("traceback replay mismatch")
    raw = args.cohort.read_bytes()
    tips = sum(line.startswith(b">") for line in raw.splitlines())
    input_hash = hashlib.sha256(raw).hexdigest()
    _, oracle = matrix(args.workspace, args.cohort, "edlib")
    out_hash = hashlib.sha256(oracle).hexdigest()
    with args.output.open("x", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("job", "case", "trial", "version", "records", "pairs", "seconds",
                         "input_sha256", "matrix_sha256", "reference"))
        for trial in range(1, 4):
            variants = [("baseline", args.baseline), ("workspace", args.workspace)]
            for name, exe in (variants if trial % 2 else variants[::-1]):
                t, result = matrix(exe, args.cohort, "nvpd-cuda")
                if result != oracle:
                    raise ValueError(f"{name} returned a different matrix")
                writer.writerow((os.environ.get("SLURM_JOB_ID", "local"), "phix174", trial, name, tips,
                                 tips*(tips-1)//2, t, input_hash, out_hash, "edlib"))
                stream.flush()


if __name__ == "__main__":
    main()
