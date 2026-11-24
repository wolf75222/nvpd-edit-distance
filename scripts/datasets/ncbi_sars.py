#!/usr/bin/env python3
"""Download the 64 versioned accessions from the published surveillance panel."""
import argparse
import csv
import hashlib
import subprocess
import tempfile
import textwrap
from collections import Counter
from pathlib import Path

from ncbi_metadata import LINEAGES, REFERENCE
from ncbi_package import read_package


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", default="datasets")
    parser.add_argument("--fasta", type=Path, default=Path("data/sars/cohort.fa"))
    args = parser.parse_args()
    base = Path(__file__).resolve().parents[2] / "results/ncbi"
    manifest = read_tsv(base / "sars_metadata.tsv")
    content = read_tsv(base / "sars_content.tsv")
    ids = [row["id"] for row in manifest]
    group_sizes = Counter(row["group"] for row in manifest)
    if (len(set(ids)) != 64 or group_sizes != {"reference": 1, **dict.fromkeys(LINEAGES, 9)}
            or [row["id"] for row in manifest if row["group"] == "reference"] != [REFERENCE]
            or ids != [row["id"] for row in content]):
        raise ValueError("invalid surveillance manifest")

    sequences = {}
    with tempfile.TemporaryDirectory() as temp:
        package = Path(temp) / "genomes.zip"
        for group in group_sizes:
            selected = [row for row in manifest if row["group"] == group]
            accessions = [row["id"] for row in selected]
            cmd = [args.datasets, "download", "virus", "genome", "accession", *accessions,
                   "--complete-only", "--include", "genome", "--no-progressbar", "--filename", str(package)]
            if group != "reference":
                cmd.extend(["--lineage", group])
            subprocess.run(cmd, check=True)
            metadata, records = read_package(package, 2697049)
            package.unlink()
            if set(accessions) != set(metadata) or set(accessions) != set(records):
                raise ValueError(f"incomplete NCBI package for {group}")
            for row in selected:
                accession = row["id"]
                record = metadata[accession]
                if (record["host"]["taxId"] != 9606 or record["virus"]["pangolinClassification"] != row["pango"]
                        or (group != "reference" and row["pango"] != group)
                        or record["length"] != int(row["nt"])
                        or len(records[accession]) != record["length"]):
                    raise ValueError(f"metadata mismatch for {accession}")
            sequences.update(records)

    fasta = []
    for row in content:
        accession = row["id"]
        seq = sequences[accession]
        counts = Counter(seq)
        ambiguous = sum(n for c, n in counts.items() if c not in "ACGT")
        if (counts.keys() - set("ACGTRYSWKMBDHVN") or counts["N"] != int(row["n"])
                or ambiguous != int(row["ambiguous"]) or hashlib.sha256(seq.encode()).hexdigest() != row["sha256"]):
            raise ValueError(f"sequence mismatch for {accession}")
        fasta.extend([">" + accession, *textwrap.wrap(seq, 80)])
    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    args.fasta.write_text("\n".join(fasta) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
