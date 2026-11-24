#!/usr/bin/env python3
"""Cross-check the full distance matrix and time each backend on a supplied FASTA."""
import argparse
import csv
import filecmp
import os
import subprocess
import sys
import time
from pathlib import Path

from run import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--fasta", required=True, type=Path)
    p.add_argument("--groups", type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--threads", required=True, type=int)
    p.add_argument("--verify-with", choices=("nvpd-cuda", "nvpd-omp"),
                   help="cross-check every pair with this NvPD implementation as well")
    args = p.parse_args()
    if args.threads < 1:
        p.error("threads must be positive")
    fasta_hash = digest(args.fasta)
    args.output.mkdir(parents=True, exist_ok=False)
    matrix = args.output / "matrix.tsv"
    algorithms = ["edlib", "wfa2"]
    if args.verify_with:
        algorithms.append(args.verify_with)
    with (args.output / "timings.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("job", "algorithm", "threads", "fasta_sha256", "matrix_sha256", "wall_seconds"))
        for algorithm in algorithms:
            path = matrix if algorithm == "edlib" else args.output / (algorithm + "-check.tsv")
            cmd = [str(args.binary.resolve()), "matrix", "--algorithm", algorithm, "--input",
                   str(args.fasta), "--threads", str(args.threads)]
            with path.open("wb") as target:
                start = time.monotonic()
                subprocess.run(cmd, check=True, stdout=target)
                elapsed = time.monotonic() - start
            writer.writerow((os.environ.get("SLURM_JOB_ID", "local"), algorithm, args.threads,
                             fasta_hash, digest(path), elapsed))
            stream.flush()
            if path != matrix:
                if not filecmp.cmp(matrix, path, shallow=False):
                    p.exit(2, f"matrices disagree. Inspect {matrix} and {path}.\n")
                path.unlink()
    tree_script = Path(__file__).resolve().parents[1] / "scripts" / "tree.py"
    for method in ("upgma", "nj"):
        cmd = [sys.executable, str(tree_script), str(matrix), "--method", method,
               "--out", str(args.output / (method + ".nwk"))]
        if args.groups:
            with (args.output / (method + "-groups.tsv")).open("w") as report:
                subprocess.run([*cmd, "--groups", str(args.groups)], check=True, stderr=report)
        else:
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
