"""Compute the exact optimum subject to |i-j| <= radius on the two periodic witnesses."""
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fasta_records import read_pair

def band_distance(a, b, radius):
    inf = len(a) + len(b) + 1
    prev = {j: j for j in range(min(len(b), radius) + 1)}
    for i, x in enumerate(a, 1):
        row = {0: i} if i <= radius else {}
        for j in range(max(1, i-radius), min(len(b), i+radius)+1):
            row[j] = min(prev.get(j-1, inf)+(x != b[j-1]),
                         prev.get(j, inf)+1, row.get(j-1, inf)+1)
        prev = row
    return prev.get(len(b), inf)

if __name__ == "__main__":
    inputs, output = map(Path, sys.argv[1:])
    with output.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("case", "radius", "constrained_distance", "unconstrained_distance"))
        for name in ("periodic512", "periodic256-tail"):
            a, b = read_pair(inputs / (name+".fa"))
            d = band_distance(a, b, max(len(a),len(b)))
            assert d == 128
            for r in (0, 16, 32, 48, 63, 64, 128):
                w.writerow((name, r, band_distance(a, b, r), d))
