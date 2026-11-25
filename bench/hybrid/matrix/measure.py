#!/usr/bin/env python3
"""Measure complete symmetric score matrices on a single GH200 allocation."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import random
import subprocess

POLICIES = ("cpu-edlib", "cpu-wfa2", "gpu", "static", "cpu-sketch", "hybrid", "hybrid-share")
COHORTS = ("near", "mixed", "distant")
FIELDS = ("cohort", "block", "position", "policy", "threads", "seconds", "sketch_seconds",
          "cpu_pairs", "gpu_pairs", "records", "pairs", "matrix_sha256")


def check_matrix(path, ids, lengths):
    with path.open(newline="") as stream:
        header, *rows = csv.reader(stream, delimiter="\t")
    if header[1:] != ids or [r[0] for r in rows] != ids:
        raise ValueError("Matrix labels do not match FASTA order")
    matrix = [list(map(int, r[1:])) for r in rows]
    if any(len(row) != len(ids) for row in matrix):
        raise ValueError("Incomplete score matrix")
    for i, row in enumerate(matrix):
        if row[i] != 0:
            raise ValueError("Nonzero diagonal")
        for j, score in enumerate(row[:i]):
            if score != matrix[j][i] or not abs(lengths[i]-lengths[j]) <= score <= max(lengths[i], lengths[j]):
                raise ValueError(f"Invalid score at ({i}, {j})")


def fasta_shape(path):
    ids, lengths = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            ids.append(line[1:].split()[0])
            lengths.append(0)
        elif line.strip():
            lengths[-1] += len(line.strip())
    return ids, lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("binary", "inputs", "output", "temp"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--panel", choices=("human", "surveillance", "composition"), default="human")
    args = parser.parse_args()
    cohorts = COHORTS if args.panel == "human" else ("surveillance",)
    panels = {
        "human": (POLICIES, 7),
        "surveillance": (("cpu-edlib", "cpu-wfa2", "hybrid"), 6),
        "composition": (("cpu-edlib", "cpu-wfa2", "hybrid", "hybrid-counts"), 8),
    }
    policies, blocks = panels[args.panel]
    args.output.mkdir(parents=True, exist_ok=True)
    args.temp.mkdir(parents=True, exist_ok=True)
    oracle = {}
    shapes = {}
    orders = {}
    rng = random.Random(839826)
    # Run an independent oracle for each complete matrix before the timed crossover.
    for cohort in cohorts:
        fasta = args.inputs / (cohort + ".fa")
        shapes[cohort] = fasta_shape(fasta)
        output = args.output / (cohort + ".tsv")
        subprocess.run([args.binary, "--input", fasta, "--policy", "cpu-edlib", "--threads", "1",
                        "--warmups", "0", "--repeats", "1", "--matrix", output], check=True,
                       stdout=subprocess.DEVNULL)
        check_matrix(output, *shapes[cohort])
        oracle[cohort] = hashlib.sha256(output.read_bytes()).hexdigest()
        if cohort == "surveillance":
            published = Path(__file__).resolve().parents[3] / "results/phylogeny/sars-edit/matrix.tsv"
            if output.read_bytes() != published.read_bytes():
                raise ValueError("The surveillance matrix differs from the published reference")
        orders[cohort] = rng.sample(policies, len(policies))

    with (args.output / "trials.csv").open("x", newline="") as times, \
         (args.output / "warmups.csv").open("x", newline="") as warm:
        writers = {phase: csv.DictWriter(f, FIELDS) for phase, f in (("measure", times), ("warmup", warm))}
        for writer in writers.values():
            writer.writeheader()
        for k in range(blocks):
            for ic in range(len(cohorts)):
                cohort = cohorts[(ic + k) % len(cohorts)]
                count = len(shapes[cohort][0])
                for pos in range(len(policies)):
                    policy = orders[cohort][(pos + k) % len(policies)]
                    matrix = args.temp / "current.tsv"
                    command = [args.binary, "--input", args.inputs / (cohort + ".fa"), "--policy", policy,
                               "--threads", str(args.threads), "--warmups", "1", "--repeats", "1", "--matrix", matrix]
                    if k == 0 and policy == ("hybrid-counts" if args.panel == "composition" else "hybrid"):
                        command += ["--decisions", args.output / (cohort + "-decisions.tsv")]
                    data = subprocess.check_output(command, text=True)
                    digest = hashlib.sha256(matrix.read_bytes()).hexdigest()
                    if digest != oracle[cohort]:
                        raise ValueError(f"Scores differ for {policy} on {cohort}")
                    rows = list(csv.DictReader(io.StringIO(data)))
                    if [r["phase"] for r in rows] != ["warmup", "measure"]:
                        raise ValueError("Expected one warm-up followed by one timed matrix")
                    for r in rows:
                        pairs = count * (count - 1) // 2
                        if (int(r["sequences"]) != count or int(r["pairs"]) != pairs
                                or int(r["cpu_pairs"]) + int(r["gpu_pairs"]) != pairs
                                or r["policy"] != policy or int(r["threads"]) != args.threads
                                or not 0 <= float(r["sketch_seconds"]) <= float(r["seconds"])
                                or not float(r["seconds"]) > 0):
                            raise ValueError(f"Invalid counters: {r}")
                        writers[r["phase"]].writerow(dict(cohort=cohort, block=k+1, position=pos+1,
                            policy=policy, threads=args.threads, seconds=r["seconds"], sketch_seconds=r["sketch_seconds"],
                            cpu_pairs=r["cpu_pairs"], gpu_pairs=r["gpu_pairs"], records=count, pairs=pairs, matrix_sha256=digest))
                    times.flush()
                    warm.flush()
                print(f"block {k+1}/{blocks} {cohort}", flush=True)


if __name__ == "__main__":
    main()
