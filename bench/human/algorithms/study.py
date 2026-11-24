#!/usr/bin/env python3
"""Paired single-core global alignments on two assemblies. All positions are zero-based."""
import argparse
import csv
import hashlib
import itertools
import json
import math
import random
import re
import resource
import statistics
import subprocess
from pathlib import Path


def read_pair(path):
    with path.open("rb") as stream:
        pair = [b"".join(line.strip() for line in group)
                for header, group in itertools.groupby(stream, lambda line: line.startswith(b">"))
                if not header]
    assert len(pair) == 2
    return pair


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def cases(root, tmp):
    with (root / "results/human/manifest.tsv").open() as stream:
        manifest = {row["case"]: row for row in csv.DictReader(stream, delimiter="\t")}
    result = []
    for case in ("chr21-2m", "chr21-4m"):
        path = root / "data/human" / (case + ".fa")
        d = manifest[case]
        assert digest(path) == d["sha256"]
        result.append((path, {k: d[k] for k in ("case", "accession_a", "start_a", "end_a",
                                                  "accession_b", "start_b", "end_b", "m", "n", "sha256")}))
    a, b = read_pair(result[-1][0])
    parent = result[-1][1]
    # Intervals read from hg38ToHs1 UCSC chain 23. Verify every 41-nt boundary
    # context is identical and appears once on each side of the original 4-Mb pair.
    for name, a0, a1, b0, b1 in (
        ("locus-34m5", 34_500_000, 34_750_000, 32_882_567, 33_132_584),
        ("locus-37m0", 37_000_000, 37_250_000, 35_382_212, 35_632_209),
    ):
        segments = []
        for seq, x, y, base in ((a, a0, a1, int(parent["start_a"])),
                                (b, b0, b1, int(parent["start_b"]))):
            x, y = x - base, y - base
            for coord in (x, y):
                flank = seq[coord-20:coord+21]
                assert len(flank) == 41 and a.count(flank) == b.count(flank) == 1
            segments.append((seq[x:y], seq[x-20:x+21], seq[y-20:y+21]))
        assert segments[0][1:] == segments[1][1:]
        path = tmp / (name + ".fa")
        with path.open("wb") as stream:
            for label, (seq, _, _) in zip(("A", "B"), segments):
                stream.write(f">{name}_{label}\n".encode())
                for i in range(0, len(seq), 80):
                    stream.write(seq[i:i+80] + b"\n")
        result.append((path, dict(case=name, accession_a=parent["accession_a"], start_a=a0, end_a=a1,
                                  accession_b=parent["accession_b"], start_b=b0, end_b=b1,
                                  m=a1-a0, n=b1-b0, sha256=digest(path))))
    return result


def replay(a, b, cigar, score):
    pos = i = j = cost = 0
    for match in re.finditer(rb"([1-9][0-9]*)([=XID])", cigar):
        run, op = int(match[1]), match[2]
        assert pos == match.start()
        pos = match.end()
        if op in (b"=", b"X"):
            left, right = a[i:i+run], b[j:j+run]
            assert len(left) == len(right) == run
            assert left == right if op == b"=" else all(x != y for x, y in zip(left, right))
        i += run if op != b"I" else 0
        j += run if op != b"D" else 0
        cost += run if op != b"=" else 0
    assert (pos, i, j, cost) == (len(cigar), len(a), len(b), score)


def cap():
    resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))


def execute(args):
    return subprocess.run(args, capture_output=True, check=True, timeout=150, preexec_fn=cap)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("build", "out", "tmp"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    inputs = cases(Path.cwd(), args.tmp)
    table(args.out / "cases.csv", [row for _, row in inputs])
    pairs, oracle = {}, {}
    for path, d in inputs:
        name = d["case"]
        pairs[name] = read_pair(path)
        ref = execute([args.build / "editdist", "pair", "--algorithm", "edlib", "--input", path])
        oracle[name] = json.loads(ref.stdout)["distance"]
    probe_input = inputs[0][0]
    counts = []
    for alg, scope in itertools.product(("wfa2", "biwfa"), ("score", "trace")):
        run = execute([args.tmp / "probe", probe_input, alg, scope])
        entries, child_bounds = map(int, run.stderr.strip().split(b","))
        assert entries == (2 if alg == "biwfa" else 0)
        assert (child_bounds > 0) == (alg == "biwfa" and scope == "trace")
        counts.append(dict(case=inputs[0][1]["case"], algorithm=alg, scope=scope,
                           calls=2, biwfa_entries=entries, recursive_bounds=child_bounds))
    table(args.out / "dispatch.csv", counts)

    rng = random.Random(260126)
    rows = []
    for block in range(1, 6):
        configs = list(itertools.product(inputs, ("biwfa", "wfa2", "edlib"), ("score", "trace")))
        rng.shuffle(configs)
        for order, ((path, d), alg, scope) in enumerate(configs, 1):
            run = execute([args.tmp / "trial", path, alg, scope])
            header, cigar = run.stdout.splitlines()
            seconds, rss, distance = header.decode().split()
            distance = int(distance)
            assert distance == oracle[d["case"]]
            if scope == "trace":
                replay(*pairs[d["case"]], cigar, distance)
            else:
                assert not cigar
            rows.append(dict(case=d["case"], algorithm=alg, scope=scope, block=block, order=order,
                             seconds=float(seconds), peak_rss_kib=int(rss), distance=distance,
                             cigar_bytes=len(cigar)))
            table(args.out / "trials.csv", rows)
    summary = []
    for key in sorted({(r["case"], r["algorithm"], r["scope"]) for r in rows}):
        group = [r for r in rows if (r["case"], r["algorithm"], r["scope"]) == key]
        record = dict(case=key[0], algorithm=key[1], scope=key[2], distance=group[0]["distance"], repeats=len(group))
        for field in ("seconds", "peak_rss_kib"):
            values = [r[field] for r in group]
            sd = statistics.stdev(values)
            record[field+"_mean"] = statistics.mean(values)
            record[field+"_sd"] = sd
            record[field+"_ci95_half"] = 2.7764451 * sd / math.sqrt(5)
        summary.append(record)
    table(args.out / "summary.csv", summary)


if __name__ == "__main__":
    main()
