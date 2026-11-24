#!/usr/bin/env python3
"""Replay the consortium's CHM13 v1.0-to-v1.1 VCF on the whole-chromosome panel."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
from pathlib import Path

PATCH_SHA256 = "bd61eca7ef07644039323f30072a463c63c7cedd24b858cc0ae582ebf61bfadc"


def read_variants(path, chromosomes):
    if hashlib.sha256(path.read_bytes()).hexdigest() != PATCH_SHA256:
        raise ValueError("unexpected v1.0-to-v1.1 release VCF")
    variants = defaultdict(list)
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alts, _, status, _, fields, sample = line.rstrip().split("\t")
            if chrom in chromosomes:
                genotype = dict(zip(fields.split(":"), sample.split(":")))["GT"]
                first = int(genotype.replace("|", "/").split("/")[0])
                alt = [ref, *alts.split(",")][first]
                if status != "PASS" or min(len(ref), len(alt)) != 1 or set(ref + alt) - set("ACGT"):
                    raise ValueError(f"unsupported release allele at {chrom}:{pos}")
                variants[chrom].append((int(pos) - 1, ref.encode(), alt.encode()))
    return variants


def check_chromosome(path, variants, m, n):
    source = bytearray()
    target = bytearray()
    with path.open("rb") as fasta:
        if fasta.readline().split() != [f">CHM13v1.0_{path.stem}".encode()]:
            raise ValueError(f"wrong first record in {path}")
        for line in fasta:
            if line.startswith(b">"):
                break
            source.extend(line.rstrip())
        if line.split() != [f">CHM13v1.1_{path.stem}".encode()] or len(source) != m:
            raise ValueError(f"wrong record or length in {path}")
        for line in fasta:
            bases = line.rstrip()
            target.extend(bases)
        if len(target) != n:
            raise ValueError(f"wrong target length in {path}")
    refseq = memoryview(source)
    patched = hashlib.sha256()
    observed = memoryview(target)
    residual_positions = []
    out = 0

    def emit(sequence):
        nonlocal out
        if out + len(sequence) > n:
            raise ValueError("patched chromosome exceeds target")
        patched.update(sequence)
        for i in range(0, len(sequence), 65536):
            block = sequence[i:i + 65536]
            ref_block = observed[out + i:out + i + len(block)]
            if block != ref_block:
                residual_positions.extend(out + i + j for j, (a, b) in enumerate(zip(block, ref_block)) if a != b)
        out += len(sequence)

    cursor = cost = near_end = end_cost = 0
    for pos, ref, alt in variants:
        end = pos + len(ref)
        if pos < cursor or refseq[pos:end] != ref:
            raise ValueError(f"VCF mismatch at {path.stem}:{pos+1}")
        emit(refseq[cursor:pos])
        emit(alt)
        # A one-letter allele aligns to one occurrence in the longer allele, if present.
        short, long = sorted((ref, alt), key=len)
        allele_cost = len(long) - 1 + int(short not in long)
        cost += allele_cost
        if end <= 1000 or pos >= m - 1000:
            near_end += 1
            end_cost += allele_cost
        cursor = end
    emit(refseq[cursor:])
    if out != n:
        raise ValueError(f"patched length differs from release on {path.stem}")
    return dict(records=len(variants), path_cost=cost, terminal_records=near_end,
                terminal_cost=end_cost, residual_substitutions=len(residual_positions),
                residual_target_positions_0=" ".join(map(str, residual_positions)),
                upper_bound=cost + len(residual_positions), patched_sha256=patched.hexdigest(),
                target_sha256=hashlib.sha256(target).hexdigest())


def verify_consensus(path, rows):
    expected = {row["case"]: row["patched_sha256"] for row in rows}
    observed = {}
    digest = None
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                name = line[1:].split()[0].decode()
                digest = None
                if name in expected:
                    if name in observed:
                        raise ValueError(f"duplicate consensus sequence: {name}")
                    digest = hashlib.sha256()
                    observed[name] = digest
            elif digest is not None:
                digest.update(line.rstrip().upper())
    if {name: digest.hexdigest() for name, digest in observed.items()} != expected:
        raise ValueError("external consensus differs from VCF replay")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", type=Path)
    p.add_argument("patch", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--consensus", type=Path, help="also verify a FASTA rebuilt by bcftools consensus -H 1")
    args = p.parse_args()
    with (args.panel / "manifest.tsv").open() as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    variants = read_variants(args.patch, {row["case"] for row in manifest})
    rows = []
    for item in manifest:
        chrom = item["case"]
        row = check_chromosome(args.panel / f"{chrom}.fa", variants[chrom], int(item["m"]), int(item["n"]))
        rows.append(dict(case=chrom, **row))
    if args.consensus:
        verify_consensus(args.consensus, rows)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
