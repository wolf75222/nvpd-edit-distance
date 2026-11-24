#!/usr/bin/env python3
"""Record a 64-sample SARS-CoV-2 panel using NCBI Datasets metadata reports."""
import argparse
import csv
import json
import subprocess
from pathlib import Path

LINEAGES = ("A", "B.1", "B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "B.1.1.529")
REFERENCE = "NC_045512.2"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", default="datasets", help="NCBI Datasets executable")
    parser.add_argument("--out", type=Path, default=Path("results/ncbi/sars_metadata.tsv"))
    args = parser.parse_args()

    def summary(*query):
        cmd = [args.datasets, "summary", "virus", "genome", *query,
               "--as-json-lines", "--complete-only"]
        lines = subprocess.check_output(cmd, text=True).splitlines()
        return [json.loads(line) for line in lines if line.strip()]

    def fields(record, group):
        if record["completeness"] != "COMPLETE" or record["virus"]["tax_id"] != 2697049:
            raise ValueError(f"unexpected metadata for {record['accession']}")
        return dict(id=record["accession"], group=group,
                    pango=record["virus"].get("pangolin_classification", ""),
                    nt=record["length"], release_date=record["release_date"][:10])

    [ref] = summary("accession", REFERENCE)
    if ref["accession"] != REFERENCE:
        raise ValueError("reference accession mismatch")
    panel = [fields(ref, "reference")]
    used = {REFERENCE}
    for lineage in LINEAGES:
        pool = summary("taxon", "2697049", "--host", "human", "--lineage", lineage,
                       "--limit", "100")
        pool.sort(key=lambda r: r["accession"])
        count = 0
        for record in pool:
            if record["accession"] in used or record["virus"].get("pangolin_classification") != lineage:
                continue
            panel.append(fields(record, lineage))
            used.add(record["accession"])
            count += 1
            if count == 9:
                break
        if count != 9:
            raise ValueError(f"only {count} distinct exact-lineage records found for {lineage}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, panel[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(panel)


if __name__ == "__main__":
    main()
