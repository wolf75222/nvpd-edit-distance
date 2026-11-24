"""Construct two orders of the same nested public chr21 windows, without editing bases."""
import argparse
import bisect
import csv
import hashlib
from pathlib import Path
import sys

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "human/algorithms"))
sys.path.insert(0, str(BENCH / "mpi/human"))
from study import read_pair
from windows import alignment_blocks, REFS

LENGTHS = (512, 2048, 4096, 12288)


def generate(assemblies, cache, output, manifest):
    a, b = read_pair(assemblies)
    for seq, key in ((a, "hg38"), (b, "hs1")):
        if hashlib.sha256(seq).hexdigest() != REFS[key][1]:
            raise ValueError(f"unexpected chromosome assembly: {key}")
    blocks = alignment_blocks(cache)
    starts = [lo for lo, _, _ in blocks]

    def anchor(x):
        pos = bisect.bisect_right(starts, x) - 1
        if pos >= 0:
            lo, hi, q = blocks[pos]
            y = q + x - lo
            if lo <= x - 20 and x + 21 <= hi and a[x-20:x+21] == b[y-20:y+21]:
                return y
        return None

    records = []
    for locus in range(16):
        origin = 19_000_000 + locus * 1_500_000
        for start in range(origin, origin + 16384, 13):
            ends = [anchor(start + n) for n in (0, *LENGTHS)]
            if None in ends:
                continue
            if b"N" not in a[start:start+max(LENGTHS)] and b"N" not in b[ends[0]:ends[-1]]:
                break
        else:
            raise ValueError(f"cannot identify all anchored ends for locus {locus}")
        for i, length in enumerate(LENGTHS,1):
            for assembly, seq, key, lo, hi in (("A",a,"hg38",start,start+length),
                                               ("B",b,"hs1",ends[0],ends[i])):
                fragment = seq[lo:hi]
                name = f"loc{locus+1:02}_{length:05}_{assembly}"
                records.append((dict(id=name, locus=locus+1, grch38_length=length, accession=REFS[key][0],
                                     start=lo, end=hi, length=len(fragment),
                                     sha256=hashlib.sha256(fragment).hexdigest()),fragment))
    output.mkdir(exist_ok=True, parents=True)
    manifest.parent.mkdir(exist_ok=True, parents=True)
    with manifest.open("w",newline="") as f:
        writer=csv.DictWriter(f,records[0][0].keys(),delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(r for r,_ in records)
    orders={"locus":records, "length":sorted(records,key=lambda v:(v[0]["length"],v[0]["id"]))}
    for name, ordered in orders.items():
        with (output/(name+".fa")).open("wb") as f:
            for r, bases in ordered:
                f.write(b">"+r["id"].encode()+b"\n"+bases+b"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("assemblies","cache","output","manifest"):
        parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    generate(args.assemblies,args.cache,args.output,args.manifest)
