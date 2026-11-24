#!/usr/bin/env python3
"""Prepare 64 complete Escherichia phage phiX174 records for the all-pairs run."""
import argparse
import csv
import hashlib
import subprocess
import tempfile
import textwrap
from pathlib import Path

from ncbi_package import read_package


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", default="datasets", help="NCBI Datasets executable")
    parser.add_argument("--fasta", type=Path, default=Path("data/phix174/cohort.fa"))
    parser.add_argument("--manifest", type=Path, default=Path("results/ncbi/phix174.tsv"),
                        help="select a new panel if absent, verify this panel if present")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        package = Path(tmp) / "phix174.zip"
        subprocess.run([args.datasets, "download", "virus", "genome", "taxon", "2886930",
                        "--complete-only", "--include", "genome", "--no-progressbar",
                        "--filename", str(package)], check=True)
        metadata, sequences = read_package(package, 2886930)

    locked = None
    if args.manifest.exists():
        with args.manifest.open() as handle:
            locked = list(csv.DictReader(handle, delimiter="\t"))
        accessions = [row["id"] for row in locked]
    else:
        reference = "NC_001422.1"
        accessions = [reference] + sorted(set(metadata) - {reference})[:63]
    if len(set(accessions)) != 64 or len(accessions) != 64:
        raise ValueError("expected 64 distinct records")

    rows, fasta = [], []
    for accession in accessions:
        record, sequence = metadata[accession], sequences[accession]
        if len(sequence) != record["length"] or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"invalid sequence for {accession}")
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        rows.append(dict(id=accession, nt=str(len(sequence)), sha256=digest))
        fasta.extend([">" + accession, *textwrap.wrap(sequence, 80)])
    if locked is not None and rows != locked:
        raise ValueError("download differs from the locked panel")
    if locked is None:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        with args.manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    args.fasta.write_text("\n".join(fasta) + "\n")


if __name__ == "__main__":
    main()
