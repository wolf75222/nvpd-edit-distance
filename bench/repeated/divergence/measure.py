#!/usr/bin/env python3
"""Paired measurements on the 60 synthetic inputs."""
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys

build, out = map(Path, sys.argv[1:])
with open("data/synthetic/manifest.tsv", newline="") as f:
    inputs = list(csv.DictReader(f, delimiter="\t"))
oracle = {}
for name in ("results/romeo/cpu-837264.csv", "results/romeo/gpu-837263.csv", "results/romeo/completion-1m/references.csv"):
    with open(name, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("status") == "ok" and row.get("distance", "").isdigit():
                previous = oracle.setdefault(row["path"], row["distance"])
                if previous != row["distance"]:
                    raise ValueError("Conflicting oracle scores: " + row["path"])
methods = ("nvpd-cuda", "edlib", "wfa2")
fields = ("job", "path", "sha256", "length", "rate", "rep", "block", "order", "algorithm", "warmup", "distance", "seconds", "rss_kib", "device_bytes")
with out.open("x", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for item in inputs:
        path = Path("data/synthetic") / item["path"]
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        expected = oracle[item["path"]]
        print(item["path"], flush=True)
        for block in range(1, 4):
            for order in range(3):
                method = methods[(order + block - 1) % 3]
                lines = subprocess.check_output([build / "repeat_kernel", path, method, expected], input="run\nrun\n", text=True).splitlines()
                if len(lines) != 3 or lines[0] != "ready":
                    raise ValueError((item["path"], method, lines))
                for warmup, line in enumerate(lines[1:]):
                    m, n, distance, seconds, rss, device = line.split(",")
                    if distance != expected or int(m) != int(item["length"]) or int(n) != int(item["mutant_length"]):
                        raise ValueError((item["path"], method, line))
                    writer.writerow(dict(job=os.environ["SLURM_JOB_ID"], path=item["path"], sha256=sha, length=item["length"], rate=item["rate"], rep=item["rep"], block=block, order=order, algorithm=method, warmup=1-warmup, distance=distance, seconds=seconds, rss_kib=rss, device_bytes=device))
                f.flush()
