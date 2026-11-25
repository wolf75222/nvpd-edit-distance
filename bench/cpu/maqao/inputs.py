#!/usr/bin/env python3
"""Extract chromosome 21 prefixes and score them independently of NvPD."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess


def prepare(inputs, binary, output):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for cohort in ("near", "distant"):
        records = (inputs / (cohort + ".fa")).read_text().split(">")[1:3]
        parsed = [rec.splitlines() for rec in records]
        sequences = ["".join(lines[1:]) for lines in parsed]
        if len(sequences) != 2 or min(map(len, sequences)) < 100000:
            raise ValueError(f"{cohort}: expected at least two windows of 100000 bases")
        for length in (32, 128, 512, 1024, 4096, 16384, 24001, 32768, 40003, 60007, 80009, 100000):
            pair = output / f"{cohort}_{length}.fa"
            pair.write_text("".join(f">{lines[0]}\n{seq[:length]}\n"
                                    for lines, seq in zip(parsed, sequences)))
            score = json.loads(subprocess.check_output(
                [binary, "pair", "--algorithm", "edlib", "--input", pair], text=True))["distance"]
            rows.append((cohort, length, parsed[0][0], parsed[1][0], score,
                         hashlib.sha256(pair.read_bytes()).hexdigest()))
    with (output / "inputs.csv").open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("cohort", "length", "first_record", "second_record", "score", "fasta_sha256"))
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "edlib_binary", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    prepare(args.inputs, args.edlib_binary, args.output)
