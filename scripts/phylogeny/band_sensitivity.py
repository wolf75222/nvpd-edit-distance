#!/usr/bin/env python3
"""Compare fixed-band edit trees to exact trees. Save aggregate counts only."""
import argparse
import ctypes as ct
import itertools
import random
import statistics
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from phylo_fit import path_matrix, write_table
from tree import build_tree, read_matrix


def read_fasta(stream):
    records, name, bases = {}, None, bytearray()
    for line in stream:
        line = line.strip()
        if not line:
            continue
        if line.startswith(b">"):
            if name is not None:
                records[name] = bytes(bases)
            fields = line[1:].split()
            if not fields:
                raise ValueError("empty FASTA identifier")
            name = fields[0].decode("ascii")
            if name in records:
                raise ValueError("duplicate FASTA identifier")
            bases.clear()
        elif name is None:
            raise ValueError("missing FASTA header")
        else:
            bases.extend(b"".join(line.upper().split()))
    if name is not None:
        records[name] = bytes(bases)
    return records


def synthetic(seed):
    """Balanced descent with substitutions and compensating string shifts."""
    rng = random.Random(seed)
    alphabet = b"0123"
    level = [bytes(rng.choices(alphabet, k=384))]
    for _ in range(5):
        children = []
        for parent in level:
            for _ in range(2):
                sequence = bytearray(parent)
                for pos in rng.sample(range(len(parent)), 4):
                    sequence[pos] = rng.choice([c for c in alphabet if c != sequence[pos]])
                width = rng.randint(1, 24)
                offset = rng.randrange(len(sequence) - width + 1)
                del sequence[offset:offset + width]
                offset = rng.randrange(len(sequence) + 1)
                sequence[offset:offset] = bytes(rng.choices(alphabet, k=width))
                children.append(bytes(sequence))
        level = children
    labeled = [(f"S{i:02}", sequence) for i, sequence in enumerate(level)]
    rng.shuffle(labeled)
    return dict(labeled)


def native_functions(path):
    lib = ct.CDLL(str(path.resolve()))
    arguments = [ct.c_char_p, ct.c_int, ct.c_char_p, ct.c_int]
    lib.optimal_witness.argtypes = arguments + [ct.POINTER(ct.c_int)]
    lib.fixed_band.argtypes = arguments + [ct.c_int]
    lib.optimal_witness.restype = lib.fixed_band.restype = ct.c_int
    return lib


def pair_costs(job):
    lib, records, reference, widths, i, j = job
    a, b = records[i], records[j]
    reach = ct.c_int()
    args = a, len(a), b, len(b)
    exact = lib.optimal_witness(*args, ct.byref(reach))
    if exact < 0 or (reference is not None and exact != reference[i][j]):
        raise ValueError(f"unrestricted oracle mismatch at row {i}, column {j}")
    endpoint = abs(len(a) - len(b))
    monotone_upper = max(len(a), len(b))
    costs = {}
    for width in sorted({max(w, endpoint) for w in widths}):
        if width >= reach.value or monotone_upper == exact:
            cost = exact
        else:
            cost = lib.fixed_band(*args, width)
        if not exact <= cost <= monotone_upper:
            raise ValueError(f"inconsistent fixed-band cost at {i}, {j}, {width}")
        costs[width] = monotone_upper = cost
    return i, j, exact, endpoint, costs


def partitions(labels, matrix, method):
    _, edges = path_matrix(build_tree(labels, matrix, method), labels)
    universe = frozenset(range(len(labels)))
    if method == "upgma":
        return {s for s, _ in edges if 1 < len(s) < len(labels)}
    return {frozenset((s, universe - s)) for s, _ in edges if 1 < len(s) < len(labels) - 1}


def run(study, native, labels, records, reference, widths, workers):
    jobs = ((native, records, reference, widths, i, j)
            for i, j in itertools.combinations(range(len(labels)), 2))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pairs = list(pool.map(pair_costs, jobs))
    modes = {"fixed": list(range(len(labels)))}
    lengths = [len(s) for s in records]
    if len(set(lengths)) > 1:
        modes["endpoint_floor"] = modes["fixed"]
        common = min(Counter(lengths).items(), key=lambda v: (-v[1], v[0]))[0]
        group = [i for i, n in enumerate(lengths) if n == common]
        if len(group) >= 4:
            modes["equal_length"] = group
    scores, trees = [], []
    for policy, members in modes.items():
        keys = [labels[i] for i in members]
        indices = {idx: pos for pos, idx in enumerate(members)}
        local = [(indices[i], indices[j], d, end, vals)
                 for i, j, d, end, vals in pairs if i in indices and j in indices]
        n = len(keys)
        matrix = [[0] * n for _ in keys]
        for i, j, d, _, _ in local:
            matrix[i][j] = matrix[j][i] = d
        truth = {m: partitions(keys, matrix, m) for m in ("upgma", "nj")}
        for width in widths:
            estimated = [[0] * n for _ in keys]
            residuals, widened = [], 0
            for i, j, d, end, vals in local:
                if policy != "endpoint_floor" and width < end:
                    continue
                widened += width < end
                estimated[i][j] = estimated[j][i] = vals[max(width, end)]
                residuals.append(estimated[i][j] - d)
            missing = len(local) - len(residuals)
            base = [study, policy, width, n]
            scores.append(base + [len(local), len(residuals), missing, widened,
                                  sum(r > 0 for r in residuals),
                                  statistics.mean(residuals) if residuals else "",
                                  max(residuals, default="")])
            if missing:
                continue
            for method, original in truth.items():
                other = partitions(keys, estimated, method)
                distance = len(original ^ other)
                denominator = len(original) + len(other)
                trees.append(base + [method, len(original), len(other), len(original & other),
                                     distance, distance / denominator])
    return scores, trees


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native", type=Path, required=True)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--fasta", help="FASTA file or - for standard input")
    source.add_argument("--synthetic", type=int, help="seed for a 32-string neutral panel")
    p.add_argument("--replicates", type=int, default=1, help="consecutive synthetic seeds")
    p.add_argument("--reference", type=Path, help="independently published square TSV")
    p.add_argument("--study", required=True, help="label written to aggregate tables")
    p.add_argument("--widths", type=int, nargs="+", required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    if args.workers < 1 or min(args.widths) < 0:
        p.error("widths must be nonnegative and workers must be positive")
    if args.replicates < 1 or (args.replicates > 1 and (args.fasta or args.reference)):
        p.error("multiple replicates require synthetic panels and no external reference matrix")
    if args.fasta == "-":
        sequences = read_fasta(sys.stdin.buffer)
    elif args.fasta:
        with open(args.fasta, "rb") as stream:
            sequences = read_fasta(stream)
    scores, trees = [], []
    native = native_functions(args.native)
    for replicate in range(args.replicates):
        study = args.study
        if args.synthetic is not None:
            seed = args.synthetic + replicate
            study += f"-{seed}"
            sequences = synthetic(seed)
        labels, ref = read_matrix(args.reference) if args.reference else (list(sequences), None)
        if set(labels) != set(sequences) or len(labels) < 4:
            p.error("need at least four unique labels, identical to the reference when supplied")
        if not all(0 < len(s) < 2**29 for s in sequences.values()) or max(args.widths) >= 2**29:
            p.error("input lengths must be positive and lengths and widths must fit the native indices")
        part_scores, part_trees = run(study, native, labels, [sequences[k] for k in labels],
                                      ref, sorted(set(args.widths)), args.workers)
        scores.extend(part_scores)
        trees.extend(part_trees)
    args.out.mkdir(parents=True, exist_ok=True)
    prefix = "study,policy,width,taxa,"
    write_table(args.out/"scores.csv", prefix+"pairs,feasible,unreachable,widened,overestimated,mean_excess,max_excess", scores)
    write_table(args.out/"trees.csv", prefix+"method,reference_parts,band_parts,shared_parts,rf,normalized_rf", trees)


if __name__ == "__main__":
    main()
