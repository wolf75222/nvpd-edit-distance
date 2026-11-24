"""Sample disjoint loci from the second chromosome record of the CHM13 v1.0/v1.1 panel."""
import argparse
import csv
import hashlib
import random
from pathlib import Path


def generate(source, folder):
    record = b"CHM13v1.1_chr1"
    bases = bytearray()
    active = False
    with source.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                active = line[1:].split()[0] == record
            elif active:
                bases.extend(line.strip().upper())
    if len(bases) != 248_387_328:
        raise ValueError("expected the complete CHM13 v1.1 chromosome 1")
    sizes = [512, 1024, 2048, 3072] * 192
    rng = random.Random(8137)
    rng.shuffle(sizes)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "chr1.fa").open("wb") as fasta, (folder / "manifest.tsv").open("w") as out:
        manifest = csv.writer(out, delimiter="\t", lineterminator="\n")
        manifest.writerow(("id", "record", "start0", "end0", "length", "ambiguous", "sha256"))
        for i, length in enumerate(sizes):
            start = rng.randrange(len(bases)*i//len(sizes), len(bases)*(i+1)//len(sizes)-length+1)
            sequence = bases[start:start+length]
            label = f"locus{i+1:04d}"
            manifest.writerow((label, record.decode(), start, start+length, length,
                               length-sum(sequence.count(c) for c in b"ACGT"), hashlib.sha256(sequence).hexdigest()))
            fasta.write(b">" + label.encode() + b"\n" + sequence + b"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    generate(args.panel, args.output)
