#!/usr/bin/env python3
"""Prepare chromosome 21 windows to test the row-layout threshold on new matrix shapes."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts/datasets"))
from human_reference import REFS, chromosome


def prepare(cache, windows, binary, output):
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    reference = {REFS[ref][0]: chromosome(cache, ref) for ref in ("hg38", "hs1")}
    with windows.open(newline="") as f:
        positions = list(csv.DictReader(f, delimiter="\t"))
    shapes = [(n, n, False) for n in (39999, 40000, 40001, 49152, 50003, 65519, 65521, 65535, 65536, 65537, 65553, 95003, 98304, 131072)]
    shapes += [(45077, 90019, False), (45077, 90019, True), (12345, 90019, False)]
    with (output / "inputs.csv").open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("cohort", "length", "columns", "first_record", "second_record", "score", "fasta_sha256"))
        for cohort in ("near", "distant"):
            source = [r for r in positions if r["cohort"] == cohort][:2]
            for rows, columns, transpose in shapes:
                parts = []
                for z, n in zip(source, (rows, columns)):
                    start = int(z["start_0"])
                    accession = z["accession"]
                    seq = reference[accession][start:start+n]
                    if len(seq) != n or seq.translate(None, b"ACGT"):
                        raise ValueError(f"Unexpected sequence {accession} at {start} length {n}")
                    parts.append((f"{accession}_{start}_{start+n}", seq))
                label = cohort
                if rows != columns:
                    label += "-rectangle"
                if transpose:
                    label += "-transpose"
                    parts.reverse()
                height, width = (len(seq) for _, seq in parts)
                path = output / f"{label}_{height}.fa"
                path.write_bytes(b"".join(b">" + name.encode() + b"\n" + seq + b"\n" for name, seq in parts))
                score = json.loads(subprocess.check_output(
                    [binary, "pair", "--algorithm", "edlib", "--input", path], text=True))["distance"]
                writer.writerow((label, height, width, parts[0][0], parts[1][0], score,
                                 hashlib.sha256(path.read_bytes()).hexdigest()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("cache", "window_manifest", "edlib_binary", "output"):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    prepare(args.cache, args.window_manifest, args.edlib_binary, args.output)
