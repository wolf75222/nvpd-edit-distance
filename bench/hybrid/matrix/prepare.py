#!/usr/bin/env python3
"""Extract three complete matrices of human chromosome 21 windows."""

import argparse
import csv
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/datasets"))
from human_reference import REFS, alignment_blocks, chromosome


LENGTH = 100_000
# Positions refer to the forward GRCh38 chromosome, with zero-based coordinates.
DISTANT = (22_500_000, 24_000_000, 26_500_000, 28_000_000,
           30_000_000, 38_000_000, 40_000_000, 42_000_000)
COHORTS = {
    "near": [(32_000_000, ref, offset) for offset in (0, 128, 256, 384)
             for ref in ("hg38", "hs1")],
    "mixed": [(anchor, ref, offset) for offset in (0, 256)
              for ref in ("hg38", "hs1") for anchor in (34_000_000, 36_000_000)],
    "distant": [(anchor, ("hg38", "hs1")[i % 2], 0) for i, anchor in enumerate(DISTANT)],
}


def prepare(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    references = {ref: chromosome(cache, ref) for ref in ("hg38", "hs1")}
    blocks = alignment_blocks(cache)
    rows = []
    for cohort, windows in COHORTS.items():
        with (out / (cohort + ".fa")).open("wb") as fasta:
            for anchor, ref, offset in windows:
                start = anchor
                if ref == "hs1":
                    hits = [q + anchor - x for x, end, q in blocks if x <= anchor < end]
                    if len(hits) != 1:
                        raise ValueError(f"anchor {anchor} lacks a unique block on chain 23")
                    start = hits[0]
                start += offset
                end = start + LENGTH
                sequence = references[ref][start:end]
                if len(sequence) != LENGTH or sequence.translate(None, b"ACGT"):
                    raise ValueError(f"invalid window {ref}:{start}-{end}")
                accession = REFS[ref][0]
                name = f"{accession}_{start}_{end}"
                fasta.write(f">{name}\n".encode("ascii"))
                for pos in range(0, LENGTH, 80):
                    fasta.write(sequence[pos:pos + 80] + b"\n")
                rows.append((cohort, name, anchor, offset, accession, start, end,
                             sequence.count(b"G") + sequence.count(b"C"),
                             hashlib.sha256(sequence).hexdigest()))

    with (out / "manifest.tsv").open("w", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("cohort", "record", "anchor_hg38_0", "offset", "accession",
                         "start_0", "end_exclusive", "gc_bases", "sequence_sha256"))
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=ROOT / "build/hybrid-matrix/cache",
                        help="cache for the pinned NCBI references and UCSC chain")
    parser.add_argument("--out", type=Path, default=ROOT / "build/hybrid-matrix/inputs")
    args = parser.parse_args()
    prepare(args.cache, args.out)
