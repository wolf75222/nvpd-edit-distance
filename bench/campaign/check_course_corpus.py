#!/usr/bin/env python3
"""Regenerate the course FASTAs and cross-check the archived ROMEO measurements."""
import argparse
import csv
import hashlib
import subprocess
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[2]
LENGTHS = (5000, 20000, 100000, 1000000)
RATES = (Decimal(".02"), Decimal(".05"), Decimal(".10"))
EDITS = ("substitutions", "insertions", "deletions")
METHODS = ("seq", "diagonal", "nvpd-omp", "nvpd-cuda", "edlib", "wfa2")


def read_corpus(directory):
    manifest = directory / "manifest.tsv"
    if manifest.read_bytes() != (ROOT / "data/synthetic/manifest.tsv").read_bytes():
        raise ValueError("manifest differs from the published seed-2026 corpus")
    with manifest.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    expected = set(product(LENGTHS, RATES, range(1, 6)))
    observed = {(int(r["length"]), Decimal(r["rate"]), int(r["rep"])) for r in rows}
    if observed != expected or len(rows) != len(expected):
        raise ValueError("expected five distinct pairs at each of the twelve settings")
    hashes, ids, bases, edits = {}, set(), Counter(), Counter()
    for row in rows:
        path = directory / row["path"]
        raw = path.read_bytes()
        headers = [f">{path.stem}_{tag}".encode() for tag in ("ref", "mut")]
        blocks = raw.split(b">")
        if len(blocks) != 3 or blocks[0]:
            raise ValueError(f"expected two FASTA records in {path}")
        sequences = []
        for block, header in zip(blocks[1:], headers):
            name, *lines = block.splitlines()
            sequence = b"".join(lines)
            if b">" + name != header or name in ids or set(sequence) - set(b"ACGT"):
                raise ValueError(f"invalid FASTA identifier or alphabet in {path}")
            ids.add(name)
            sequences.append(sequence)
        length = int(row["length"])
        counts = {key: int(row[key]) for key in EDITS}
        mutant_length = length + counts["insertions"] - counts["deletions"]
        if (min(counts.values()) < 0 or sum(counts.values()) != length * Decimal(row["rate"])
                or mutant_length != int(row["mutant_length"])
                or list(map(len, sequences)) != [length, mutant_length]):
            raise ValueError(f"inconsistent injected edits or sequence lengths in {path}")
        bases.update(sequences[0].decode("ascii"))
        edits.update(counts)
        hashes[hashlib.sha256(raw).hexdigest()] = row
    if len(hashes) != len(rows):
        raise ValueError("duplicate FASTA content")
    return hashes, bases, edits


def check_runs(pairs):
    campaign = ROOT / "results/romeo"
    scores, coverage = {}, defaultdict(set)
    files = sorted(campaign.glob("cpu-*.csv")) + sorted(campaign.glob("gpu-*.csv"))
    files.append(campaign / "completion-1m/trials.csv")
    for path in files:
        with path.open() as stream:
            for row in csv.DictReader(stream):
                key, algorithm = row["input_sha256"], row["algorithm"]
                sample = pairs[key]
                if (row["status"] != "ok" or row["m"] != sample["length"]
                        or row["n"] != sample["mutant_length"]
                        or any(row[k] != sample[k] for k in ("path", "length", "rate", "rep", "seed"))):
                    raise ValueError(f"inconsistent ROMEO measurement for {key} in {path}")
                score = int(row["distance"])
                if not abs(int(row["m"]) - int(row["n"])) <= score <= sum(int(sample[k]) for k in EDITS):
                    raise ValueError(f"distance outside edit-count bounds for {sample['path']} in {path}")
                if scores.setdefault(key, score) != score:
                    raise ValueError(f"distance disagreement for {sample['path']} in {path}")
                coverage[algorithm].add(key)
    scalar = {key for key, r in pairs.items() if int(r["length"]) < 1000000 or r["rep"] == "1"}
    for algorithm in METHODS:
        if coverage[algorithm] != (scalar if algorithm == "seq" else pairs.keys()):
            raise ValueError(f"incomplete or unexpected ROMEO inputs for {algorithm}")
    large_scalar = {key for key in scalar if pairs[key]["length"] == "1000000"}
    with (ROOT / "results/repeated/epyc-scalar/trials.csv").open() as stream:
        repeated = list(csv.DictReader(stream))
    if ({r["input_sha256"] for r in repeated} != large_scalar
            or any(r["algorithm"] != "seq" or int(r["distance"]) != scores[r["input_sha256"]] for r in repeated)):
        raise ValueError("the repeated scalar trials differ from the three 1 Mnt reference pairs")
    return coverage


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, help="audit an existing default synthetic corpus")
    args = ap.parse_args()
    with TemporaryDirectory() as tmp:
        directory = Path(tmp) if args.data is None else args.data
        if args.data is None:
            subprocess.run([sys.executable, str(ROOT / "scripts/datasets/synthetic.py"), "--out", tmp], check=True)
        pairs, bases, edits = read_corpus(directory)
        coverage = check_runs(pairs)
    print(f"{len(pairs)} distinct pairs, {sum(edits.values())} injected edits")
    print(", ".join(f"{base}={bases[base]}" for base in "ACGT"))
    print(", ".join(f"{kind}={edits[kind]}" for kind in EDITS))
    print("Validated pairs: " + ", ".join(f"{algorithm}={len(coverage[algorithm])}" for algorithm in METHODS))


if __name__ == "__main__":
    main()
