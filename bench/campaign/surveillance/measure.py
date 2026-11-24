#!/usr/bin/env python3
"""Time whole 64-genome matrices. Reject each run unless its TSV matches the reference."""
import argparse
import csv
import hashlib
import os
import random
import subprocess
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("target", choices=("cpu", "gpu", "release"))
    parser.add_argument("--repetitions", type=int, default=5)
    args = parser.parse_args()
    if args.repetitions < 2:
        parser.error("at least two timed blocks are required")
    reference = args.reference.read_bytes()
    digest = hashlib.sha256(reference).hexdigest()
    threads = int(os.environ["SLURM_CPUS_PER_TASK"])
    if args.target == "cpu":
        configs = [(alg, width) for alg in ("edlib", "wfa2", "nvpd-adaptive", "nvpd-split")
                   for width in (1, threads)]
        configs.append(("nvpd-omp", threads))
    elif args.target == "release":
        configs = [(alg, width) for alg in ("edlib", "wfa2", "biwfa", "nvpd-split")
                   for width in (1, threads)]
    else:
        configs = [("nvpd-cuda-local", 1), ("edlib", 1), ("wfa2", 1)]
    rng = random.Random(1001)
    with args.output.open("x", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(("job", "block", "position", "algorithm", "threads", "wall_s", "matrix_sha256"))
        for block in range(args.repetitions + 1):
            rng.shuffle(configs)
            for position, (algorithm, width) in enumerate(configs, 1):
                cmd = [str(args.binary.resolve()), "matrix", "--input", str(args.fasta),
                       "--algorithm", algorithm, "--threads", str(width)]
                start = time.perf_counter_ns()
                matrix = subprocess.run(cmd, check=True, stdout=subprocess.PIPE).stdout
                elapsed = (time.perf_counter_ns() - start) / 1e9
                if matrix != reference:
                    raise ValueError(f"matrix differs from reference for {algorithm}, block {block}")
                writer.writerow((os.environ["SLURM_JOB_ID"], block, position, algorithm,
                                 width, f"{elapsed:.9f}", digest))
                out.flush()


if __name__ == "__main__":
    main()
