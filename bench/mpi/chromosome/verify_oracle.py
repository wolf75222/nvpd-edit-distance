"""Check archived MPI scores and certificate counts against a serial-format edlib matrix."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from measure import configurations, examine


def verify(matrix, cohort, results):
    with (cohort / "manifest.tsv").open() as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    fasta = cohort / "chr1.fa"
    hashes = {name: digest for digest, name in (row.split() for row in (results / "inputs.sha256").read_text().splitlines())}
    for path in (fasta, cohort / "manifest.tsv"):
        if hashes[path.name] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError(f"{path.name} differs from the ROMEO input")
    lengths = {row["id"]: int(row["length"]) for row in manifest}
    if len(lengths) != 768:
        raise ValueError("expected 768 distinct locus identifiers")
    text = matrix.read_text()
    _, eligible, matching = examine(text, lengths)
    digest = hashlib.sha256(text.encode()).hexdigest()
    with (results / "trials.csv").open() as stream:
        trials = [r for r in csv.DictReader(stream) if int(r["records"]) == len(lengths)]
    expected_trials = {(block, nodes, solver, policy, int(bound)) for block in range(1, 6)
                       for size, nodes, solver, policy, bound in configurations() if size == len(lengths)}
    actual_trials = {(int(r["trial"]), int(r["nodes"]), r["algorithm"], r["schedule"], int(r["length_bound"]))
                     for r in trials}
    if actual_trials != expected_trials or len(trials) != len(expected_trials):
        raise ValueError("incomplete or duplicated MPI trial")
    for row in trials:
        expected = matching if int(row["length_bound"]) else 0
        if row["matrix_sha256"] != digest or int(row["certified"]) != expected:
            raise ValueError(f"edlib disagrees with MPI: {row}")
    record = dict(oracle="edlib", records=len(lengths), pairs=len(lengths)*(len(lengths)-1)//2,
                  fasta_sha256=hashes[fasta.name], matrix_sha256=digest, eligible_pairs=eligible,
                  certified_pairs=matching, verified_mpi_trials=len(trials))
    (results / "edlib-oracle.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path, help="TSV from editdist matrix --algorithm edlib")
    parser.add_argument("cohort", type=Path, help="directory with chr1.fa and manifest.tsv")
    parser.add_argument("results", type=Path, help="directory containing completed trials.csv and inputs.sha256")
    args = parser.parse_args()
    verify(args.matrix, args.cohort, args.results)
