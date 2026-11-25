"""Extract two nested sets of homologous chr21 windows from chain 23."""

import argparse
from bisect import bisect_right
import csv
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts/datasets"))
from human_reference import REFS, alignment_blocks, chromosome


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path, help="directory for the pinned NCBI sequences and UCSC chain")
    p.add_argument("out", type=Path)
    args = p.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    a, b = (chromosome(args.cache, k) for k in ("hg38", "hs1"))
    blocks = alignment_blocks(args.cache)
    starts = [x[0] for x in blocks]

    def lift(k):
        i, end, j = blocks[bisect_right(starts, k) - 1]
        if not i <= k < end:
            raise ValueError("endpoint falls outside aligned chain block")
        return j + k - i

    rows = []
    reaches = []
    for region, center in [("collinear", 37400000), ("insertion", 33474775)]:
        for m in (16384, 65536, 131072, 262144, 524288):
            i = center - m // 2
            j = lift(i)
            end = lift(i + m)
            aa, bb = a[i : i + m], b[j:end]
            if set(aa) != set(b"ACGT") or set(bb) != set(b"ACGT"):
                raise ValueError("expected four canonical bases in both chromosome windows")
            name = f"{region}-{m}"
            fa = b">GRCh38\n" + aa + b"\n>CHM13\n" + bb + b"\n"
            (args.out / (name + ".fa")).write_bytes(fa)
            rows.append(
                dict(
                    case=name,
                    region=region,
                    accession_a=REFS["hg38"][0],
                    start_a=i,
                    end_a=i + m,
                    accession_b=REFS["hs1"][0],
                    start_b=j,
                    end_b=end,
                    chain_id=23,
                    strand="+",
                    m=m,
                    n=len(bb),
                    input_sha256=hashlib.sha256(fa).hexdigest(),
                )
            )
            # Native CUDA tiles depend on the farthest previous-row column the
            # recurrence consults. Before a symbol first appears only j and j-1 are used.
            pattern, text = (aa, bb) if len(aa) <= len(bb) else (bb, aa)
            reach = 1
            for symbol in set(pattern):
                k = 0
                for j, c in enumerate(text, 1):
                    if c == symbol:
                        k = j
                    if k:
                        reach = max(reach, j - k + 1)
            reaches.append(
                dict(
                    case=name,
                    input_sha256=rows[-1]["input_sha256"],
                    maximum_match_reach=reach,
                    native_cuda_path="eight-row tile" if reach <= 64 else "one row",
                )
            )
    with (args.out / "reach.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=reaches[0], lineterminator="\n")
        w.writeheader()
        w.writerows(reaches)
    with (args.out / "inputs.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
