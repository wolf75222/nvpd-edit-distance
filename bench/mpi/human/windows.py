#!/usr/bin/env python3
"""Match disjoint human chromosome 21 windows on UCSC's hg38 to hs1 chain."""
import argparse
import bisect
import csv
import hashlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts/datasets"))
from human_reference import REFS, alignment_blocks, chromosome


def generate(cache, out, manifest):
    hg, hs = (chromosome(cache, key) for key in ("hg38", "hs1"))
    blocks = alignment_blocks(cache)
    starts = [b[0] for b in blocks]

    def anchor(x):
        idx = bisect.bisect_right(starts, x) - 1
        if idx < 0:
            return None
        lo, hi, q0 = blocks[idx]
        y = q0 + x - lo
        if lo <= x - 20 and x + 21 <= hi and hg[x-20:x+21] == hs[y-20:y+21]:
            return y
        return None

    def select(origin, length):
        for x in range(origin, origin + 8_192, 13):
            y, z = anchor(x), anchor(x + length)
            if y is None or z is None or b"N" in hg[x:x+length] or b"N" in hs[y:z]:
                continue
            return x, x+length, y, z
        return None

    # Keep the original spacing, including gaps from incomplete/ambiguous sequence.
    grid = []
    for origin in range(19_000_000, 46_000_000, 28_000):
        window = select(origin, 16_384)
        if window:
            grid.append(window)
        if len(grid) == 768:
            break
    if len(grid) != 768 or any(u[1] > v[0] or u[3] > v[2] for u,v in zip(grid, grid[1:])):
        raise ValueError("expected 768 non-overlapping homologous intervals")
    # Alternate along the chromosome so the half cohort spans the same interval.
    grid = grid[::2] + grid[1::2]
    out.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="") as table, (out / "matched-768.fa").open("wb") as large, (out / "matched-384.fa").open("wb") as small:
        wr = csv.writer(table, delimiter="\t", lineterminator="\n")
        wr.writerow(("id","accession_a","start_a","end_a","accession_b","start_b","end_b","m","n","cells","sha256_a","sha256_b"))

        def records(i, w):
            x, e, y, z = w
            a, b = hg[x:e], hs[y:z]
            wr.writerow((i,REFS["hg38"][0],x,e,REFS["hs1"][0],y,z,len(a),len(b),len(a)*len(b),hashlib.sha256(a).hexdigest(),hashlib.sha256(b).hexdigest()))
            return f">{i}_A\n".encode()+a+b"\n"+f">{i}_B\n".encode()+b+b"\n"

        for j,w in enumerate(grid):
            fasta = records(f"locus{j+1:03}",w)
            large.write(fasta)
            if j < 384:
                small.write(fasta)
        for i,origin,length in [("pair-64k",30_000_000,65_536),("pair-128k",40_000_000,131_072)]:
            w = select(origin,length)
            if w is None:
                raise ValueError(f"unable to match {i}")
            (out / (i+".fa")).write_bytes(records(i,w))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("manifest",type=Path)
    opts = parser.parse_args()
    opts.cache.mkdir(parents=True, exist_ok=True)
    generate(opts.cache, opts.output, opts.manifest)
