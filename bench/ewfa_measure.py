"""Validate raw GPU scores/paths on exact input pairs, then time the original release."""
import argparse
import csv
import hashlib
import random
import re
import subprocess
import time
from pathlib import Path

from alternatives import fasta_pair


def run(cmd):
    start = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=True)
    if "Error cuda" in p.stdout or "This is a bug" in p.stdout or "Max distance reached" in p.stdout:
        raise RuntimeError(p.stdout[:400])
    return p.stdout, time.perf_counter() - start


def replay(a, b, cigar):
    i = j = cost = 0
    tokens = list(re.finditer(r"(\d*)([MXID])", cigar))
    if "".join(tok.group() for tok in tokens) != cigar:
        return None
    for tok in tokens:
        n, op = int(tok[1] or 1), tok[2]
        left, right = i + (0 if op == "I" else n), j + (0 if op == "D" else n)
        if n <= 0 or left > len(a) or right > len(b):
            return None
        if op == "M" and a[i:left] != b[j:right]:
            return None
        if op == "X" and any(x == y for x, y in zip(a[i:left], b[j:right])):
            return None
        cost += n if op != "M" else 0
        i, j = left, right
    return cost if (i, j) == (len(a), len(b)) else None


def write(path, rows):
    with path.open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def scores(text):
    return [int(d) for d in re.findall(r"^distance (\d+)$", text, re.M)]


def main(work, out):
    with (work/"synthetic/manifest.tsv").open() as f:
        manifest = list(csv.DictReader(f, delimiter="\t"))
    rows, trials = [], []
    rng = random.Random(2026)
    for case, rate, cap in [("5k_full_2pct", "0.02", 5000*2),
                            ("1k_prefix_2pct", "0.02", 1000), ("1k_prefix_5pct", "0.05", 1000),
                            ("5k_full_low", "0.005", 5000*2)]:
        source = [work/"synthetic"/r["path"] for r in manifest if r["rate"] == rate]
        pairs = [(a[:cap], b[:cap]) for a, b in map(fasta_pair, source)]
        seq = work/"case.seq"
        seq.write_text("".join(f">{a}\n<{b}\n" for a,b in pairs))
        oracle = scores(run([work/"baselines","edlib",seq])[0])
        if len(oracle) != len(pairs):
            raise RuntimeError("missing oracle scores")
        selected = [(i,a,b) for i,(a,b) in enumerate(pairs) if oracle[i] < 64]
        validated, stock = {}, {}
        if selected:
            seq.write_text("".join(f">{a}\n<{b}\n" for _,a,b in selected))
            width = max(max(len(a),len(b)) for _,a,b in selected)
            opts = ["-f",seq,"-n",str(len(selected)),"-l",str(width),
                    "-b",str(len(selected)),"-t","1"]
            raw = run([work/"raw",*opts,"-p"])[0]
            for match in re.finditer(r"^raw (\d+) (\d+) (\S*)$",raw,re.M):
                i,a,b = selected[int(match[1])]
                validated[i] = (int(match[2]), replay(a,b,match[3]))
            printed = run([work/"original",*opts,"-p"])[0]
            for match in re.finditer(r"^(\d+): (\S*)$",printed,re.M):
                i,a,b = selected[int(match[1])]
                stock[i] = replay(a,b,match[2]) == oracle[i]
            batch_valid = all(validated.get(i) == (oracle[i],oracle[i]) for i,_,_ in selected)
            algos = ["edlib","nvpd-cuda"] + (["ewfa-gpu-2022"] if batch_valid else [])
            for trial in range(6):
                rng.shuffle(algos)
                for algo in algos:
                    if algo == "ewfa-gpu-2022":
                        line, wall = run([work/"original",*opts])
                        sec = float(re.search(r"Aligned executed\. Wall time ([\d.]+)s",line)[1])
                    else:
                        line, wall = run([work/"baselines",algo,seq])
                        if scores(line) != [oracle[i] for i,_,_ in selected]:
                            raise RuntimeError(f"baseline score disagreement for {case}, {algo}")
                        sec = float(re.search(r"^seconds (\S+)",line,re.M)[1])
                    if trial:
                        trials.append(dict(case=case, algorithm=algo, trial=trial, pairs=len(selected),
                                           worker_threads=1, batch_seconds=sec, process_seconds=wall))
        for i,(a,b) in enumerate(pairs):
            raw_score, raw_cost = validated.get(i,("", ""))
            rows.append(dict(case=case,source=source[i].name,source_sha256=hashlib.sha256(source[i].read_bytes()).hexdigest(),
                             pair_sha256=hashlib.sha256(f">{a}\n<{b}\n".encode()).hexdigest(),
                             n=len(a),m=len(b),edlib=oracle[i],gpu_score=raw_score,raw_cigar_cost=raw_cost,
                             stock_printed_cigar_replays=stock.get(i,""),status="outside_d_lt_64" if oracle[i]>=64 else
                             ("verified" if (raw_score,raw_cost)==(oracle[i],oracle[i]) else "validation_failed")))
    write(out/"pairs.csv",rows)
    write(out/"trials.csv",trials)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args = p.parse_args()
    main(args.work,args.out)
