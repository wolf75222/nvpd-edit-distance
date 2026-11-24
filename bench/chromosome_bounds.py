"""Bound the global edit distance using base counts and one collinear UCSC chain."""
import argparse
from collections import Counter
import csv
import gzip
from pathlib import Path

from fasta_records import read_pair


def column_cost(a, b):
    # Each of the two end-anchored alignments gives a constructive upper bound.
    n = min(len(a), len(b))
    gap = abs(len(a) - len(b))
    if not n:
        return gap
    left = sum(x != y for x, y in zip(a[:n], b[:n]))
    if gap == 0:
        return left
    right = sum(x != y for x, y in zip(a[-n:], b[-n:]))
    return min(left, right) + gap


def bounds(a, b, chain, number):
    ca, cb = Counter(a), Counter(b)
    surplus_a, surplus_b = sum((ca - cb).values()), sum((cb - ca).values())
    a, b = memoryview(a), memoryview(b)
    with gzip.open(chain, "rt") as stream:
        for line in stream:
            fields = line.split()
            if fields[:1] == ["chain"] and int(fields[-1]) == number:
                break
        else:
            raise ValueError("chain ID absent")
        if len(fields) != 13 or fields[2] != "chr21" or fields[7] != "chr21" or fields[4] != "+" or fields[9] != "+":
            raise ValueError("expected a positive chr21/chr21 chain")
        if (int(fields[3]), int(fields[8])) != (len(a), len(b)):
            raise ValueError("incompatible chain lengths")
        i, j = int(fields[5]), int(fields[10])
        end_a, end_b = int(fields[6]), int(fields[11])
        if not (0 <= i < end_a <= len(a) and 0 <= j < end_b <= len(b)):
            raise ValueError("chain interval lies outside the FASTA")
        upper = column_cost(a[:i], b[:j])
        for line in stream:
            block = list(map(int, line.split()))
            if len(block) not in (1, 3) or block[0] <= 0 or any(v < 0 for v in block[1:]):
                raise ValueError("invalid chain block")
            size = block[0]
            da, db = block[1:] if len(block) == 3 else (0, 0)
            if i + size + da > end_a or j + size + db > end_b:
                raise ValueError("chain block extends past endpoint")
            upper += column_cost(a[i:i+size], b[j:j+size])
            i += size
            j += size
            if len(block) == 1:
                break
            upper += column_cost(a[i:i+da], b[j:j+db])
            i += da
            j += db
        else:
            raise ValueError("chain lacks a terminal block")
        if (i, j) != (end_a, end_b):
            raise ValueError("inconsistent block endpoints")
        upper += column_cost(a[i:], b[j:])
    return dict(chain=number, m=len(a), n=len(b), length_gap=abs(len(a)-len(b)),
                surplus_a=surplus_a, surplus_b=surplus_b, lower=max(surplus_a,surplus_b), upper=upper)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("chain", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    row = bounds(*read_pair(args.fasta), args.chain, 23)
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, row.keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
