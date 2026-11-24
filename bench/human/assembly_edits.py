"""List nonmatching runs of an optimal global path between two FASTA records."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from fasta_records import read_pair


def edit_runs(a, b, cigar, distance):
    """Replay from A to B; output zero-based, half-open intervals in both assemblies."""
    i = j = cost = pos = 0
    rows = []
    for match in re.finditer(r"([1-9][0-9]*)([=XID])", cigar):
        if match.start() != pos:
            raise ValueError("Invalid CIGAR")
        length, op = int(match[1]), match[2]
        end_a = i + (length if op != "I" else 0)
        end_b = j + (length if op != "D" else 0)
        if end_a > len(a) or end_b > len(b):
            raise ValueError("CIGAR extends beyond the FASTA input")
        left, right = a[i:end_a], b[j:end_b]
        if op == "=" and left != right:
            raise ValueError("Different bases inside a matching run")
        if op == "X" and any(x == y for x, y in zip(left, right)):
            raise ValueError("Equal bases inside a substitution run")
        if op != "=":
            rows.append((i, end_a, j, end_b, op, length))
            cost += length
        i, j, pos = end_a, end_b, match.end()
    if (i, j, pos, cost) != (len(a), len(b), len(cigar), distance):
        raise ValueError("CIGAR span or cost differs from alignment result")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    cmd = [args.binary.resolve(), "pair", "--algorithm", "edlib", "--trace", "--input", args.fasta]
    result = json.loads(subprocess.check_output(cmd))
    a, b = read_pair(args.fasta)
    rows = edit_runs(memoryview(a), memoryview(b), result["cigar"], result["distance"])
    with args.output.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("start_a", "end_a", "start_b", "end_b", "operation", "length"))
        writer.writerows(rows)


if __name__ == "__main__":
    main()
