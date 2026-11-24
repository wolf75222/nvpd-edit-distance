#!/usr/bin/env python3
"""Generate reproducible pairs of synthetic DNA sequences and a TSV manifest."""

import argparse
import csv
import hashlib
import os
import random
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path
from tempfile import NamedTemporaryFile


BASES = b"ACGT"
OTHER_BASES = {base: bytes(b for b in BASES if b != base) for base in BASES}
FIELDS = ("path", "length", "rate", "rep", "seed", "mutant_length",
          "substitutions", "insertions", "deletions")


def pair_seed(master, length, rate, rep):
    key = f"{master}/{length}/{rate}/{rep}".encode("ascii")
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "big")


def make_pair(length, edits, seed):
    rng = random.Random(seed)
    reference = bytes(rng.choices(BASES, k=length))
    mutant = bytearray()
    counts = [0, 0, 0]
    start = 0
    for pos in sorted(rng.sample(range(length), edits)):
        mutant.extend(reference[start:pos])
        op = rng.randrange(3)
        counts[op] += 1
        if op == 0:
            mutant.append(rng.choice(OTHER_BASES[reference[pos]]))
        elif op == 1:
            mutant.append(reference[pos])
            mutant.append(rng.choice(BASES))
        # A deletion consumes the reference base without emitting a base.
        start = pos + 1
    mutant.extend(reference[start:])
    return reference, mutant, counts


@contextmanager
def atomic_file(path, mode, **options):
    temp = NamedTemporaryFile(mode, dir=path.parent, delete=False, **options)
    try:
        with temp:
            yield temp
        os.replace(temp.name, path)
    finally:
        Path(temp.name).unlink(missing_ok=True)


def write_pair(path, pair_id, reference, mutant):
    with atomic_file(path, "wb") as out:
        for kind, seq in (("ref", reference), ("mut", mutant)):
            out.write(f">{pair_id}_{kind}\n".encode("ascii"))
            for start in range(0, len(seq), 80):
                out.write(seq[start:start + 80] + b"\n")


def generate(out, lengths, rates, pairs, master_seed):
    pair_dir = out / "pairs"
    pair_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for length in lengths:
        for rate in rates:
            rate_text = format(rate.normalize(), "f")
            edits = int((length * rate).to_integral_value(rounding=ROUND_HALF_UP))
            for rep in range(1, pairs + 1):
                seed = pair_seed(master_seed, length, rate_text, rep)
                reference, mutant, counts = make_pair(length, edits, seed)
                pair_id = f"n{length}_r{rate_text}_p{rep}_{seed:016x}"
                path = pair_dir / (pair_id + ".fa")
                write_pair(path, pair_id, reference, mutant)
                rows.append((path.relative_to(out), length, rate_text, rep, seed,
                             len(mutant), *counts))
    with atomic_file(out / "manifest.tsv", "w", encoding="utf-8", newline="") as manifest:
        writer = csv.writer(manifest, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        writer.writerows(rows)


def rate_arg(text):
    try:
        rate = Decimal(text)
    except InvalidOperation:
        raise argparse.ArgumentTypeError("expected a decimal error fraction") from None
    if not rate.is_finite() or not 0 <= rate <= 1:
        raise argparse.ArgumentTypeError("rate must be between 0 and 1")
    return rate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("data/synthetic"))
    parser.add_argument("--lengths", type=int, nargs="+", default=[5000, 20000, 100000, 1000000])
    parser.add_argument("--rates", type=rate_arg, nargs="+", default=list(map(Decimal, [".02", ".05", ".10"])))
    parser.add_argument("--pairs", type=int, default=5, help="pairs per length and rate")
    parser.add_argument("--seed", type=int, default=2026, help="master seed")
    args = parser.parse_args()
    if args.pairs < 1 or min(args.lengths) < 1:
        parser.error("lengths and pair count must be positive")
    if len(set(args.lengths)) != len(args.lengths) or len(set(args.rates)) != len(args.rates):
        parser.error("duplicate lengths or rates")
    generate(args.out, args.lengths, args.rates, args.pairs, args.seed)


if __name__ == "__main__":
    main()
