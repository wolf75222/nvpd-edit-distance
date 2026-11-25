#!/usr/bin/env python3
"""Nested chain-23 chr21q alignments, one fixed EPYC core per window length."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import math
import os
import random
import resource
import statistics
import subprocess
from pathlib import Path

from study import digest, read_pair, replay, table

START = 14_000_100
SIZES = (8, 16, 24, 32)
SEQ_HASH = ("ff47a5bff4a148bc684e444e4c52eaceec1290cb6531b2d8c613becdebbad23c",
            "f210b99cc28cb8bf66a9aca6d0cd6ceaab1ef8a78ff135876f25644147e27481")
CHAIN_HASH = "926d4583eda8f90c88434293859b06906573f2469630e898b8b475ac8751327a"


def windows(tmp, chain):
    a, b = read_pair(Path("data/human/chr21-assemblies.fa"))
    assert tuple(hashlib.sha256(x).hexdigest() for x in (a, b)) == SEQ_HASH
    assert digest(chain) == CHAIN_HASH
    hits = {START + mb*1_000_000: [] for mb in (0, *SIZES)}
    active = False
    with gzip.open(chain, "rt") as stream:
        for line in stream:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                active = f[2] == "chr21" and f[4] == "+"
                left, right = int(f[5]), int(f[10])
                strand, query, cid = f[9], f[7], f[-1]
            elif active:
                size = int(f[0])
                for x in hits:
                    if left <= x < left + size:
                        hits[x].append((cid, query, strand, right+x-left))
                if len(f) == 3:
                    left += size+int(f[1])
                    right += size+int(f[2])
    points = {}
    for x, records in hits.items():
        assert len(records) == 1 and records[0][:3] == ("23", "chr21", "+"), (x, records)
        points[x] = y = records[0][3]
        context = a[x-20:x+21]
        assert context == b[y-20:y+21] and a.count(context) == b.count(context) == 1
    cases = []
    for mb in SIZES:
        end = START + mb*1_000_000
        l, r = a[START:end], b[points[START]:points[end]]
        case = f"chr21q-{mb}m"
        path = tmp / (case + ".fa")
        with path.open("wb") as stream:
            for label, seq in zip(("GRCh38", "CHM13"), (l, r)):
                stream.write(f">{label}_{case}\n".encode())
                for i in range(0, len(seq), 80):
                    stream.write(seq[i:i+80]+b"\n")
        cases.append(dict(case=case, accession_a="NC_000021.9", start_a=START, end_a=end,
                          accession_b="NC_060945.1", start_b=points[START], end_b=points[end],
                          m=len(l), n=len(r), N_a=l.count(b"N"), N_b=r.count(b"N"),
                          sha256=digest(path)))
    return cases


def process_cap():
    resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))


def execute(argv):
    return subprocess.run(argv, capture_output=True, check=True, timeout=720, preexec_fn=process_cap)


def measure(args, case, cpu):
    os.sched_setaffinity(0, {cpu})
    assert os.sched_getaffinity(0) == {cpu}
    name = case["case"]
    path = args.tmp / (name + ".fa")
    a, b = read_pair(path)
    ref = execute([args.build / "editdist", "pair", "--algorithm", "wfa2", "--input", path])
    oracle = json.loads(ref.stdout)["distance"]
    methods = [("biwfa", "trace"), ("biwfa", "score"), ("wfa2", "score")]
    if name == "chr21q-8m":
        methods += [("edlib", "score"), ("edlib", "trace")]
    rng = random.Random(260126 + case["m"])
    rows, capped = [], set()
    for block in range(1, 6):
        order = methods.copy()
        rng.shuffle(order)
        for position, (alg, scope) in enumerate(order, 1):
            row = dict(case=name, algorithm=alg, scope=scope, block=block, order=position, cpu=cpu,
                       seconds="", peak_rss_kib="", distance="", cigar_bytes="", status="ok")
            if (alg, scope) in capped:
                row["status"] = "skipped_after_cap"
            else:
                try:
                    run = execute([args.tmp / "trial", path, alg, scope])
                except subprocess.TimeoutExpired:
                    row["status"] = "timeout_720s"
                    capped.add((alg, scope))
                else:
                    header, cigar = run.stdout.splitlines()
                    seconds, rss, distance = header.decode().split()
                    assert int(distance) == oracle
                    if scope == "trace":
                        replay(a, b, cigar, oracle)
                    else:
                        assert not cigar
                    row.update(seconds=float(seconds), peak_rss_kib=int(rss), distance=oracle,
                               cigar_bytes=len(cigar))
            rows.append(row)
            table(args.out / (name + ".partial.csv"), rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("build", "out", "tmp", "chain"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    cases = windows(args.tmp, args.chain)
    cpus = sorted(os.sched_getaffinity(0))
    assert len(cpus) == len(cases) == 4, cpus
    table(args.out / "cases.csv", cases)
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as workers:
        tasks = [workers.submit(measure, args, case, cpu) for case, cpu in zip(cases, cpus)]
        rows = [r for task in tasks for r in task.result()]
    table(args.out / "trials.csv", rows)
    summary = []
    for case, alg, scope in sorted({(r["case"], r["algorithm"], r["scope"]) for r in rows}):
        group = [r for r in rows if (r["case"], r["algorithm"], r["scope"]) == (case, alg, scope)
                 and r["status"] == "ok"]
        count = len(group)
        record = dict(case=case, algorithm=alg, scope=scope, distance=group[0]["distance"] if count else "",
                      repeats=count, incomplete=5-count)
        for field in ("seconds", "peak_rss_kib"):
            values = [r[field] for r in group]
            record[field+"_mean"] = statistics.mean(values) if count else ""
            record[field+"_sd"] = statistics.stdev(values) if count > 1 else ""
            record[field+"_ci95_half"] = 2.7764451*statistics.stdev(values)/math.sqrt(5) if count == 5 else ""
        summary.append(record)
    table(args.out / "summary.csv", summary)
    for case in cases:
        (args.out / (case["case"] + ".partial.csv")).unlink()


if __name__ == "__main__":
    main()
