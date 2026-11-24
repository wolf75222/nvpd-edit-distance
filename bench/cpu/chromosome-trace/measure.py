"""Randomized complete blocks. Each new process warms its solver before timing."""
import csv
import os
import random
import subprocess
import sys
from pathlib import Path

binary, panel, divergent, output = map(Path, sys.argv[1:])
inputs = {c: panel / f"{c}.fa" for c in ("chr1", "chr18", "chr20")}
inputs["chr21-2m"] = divergent
affinity = sorted(os.sched_getaffinity(0))
if len(affinity) < 8:
    raise ValueError(f"need eight physical cores, got {affinity}")
configs = [(case, mode, method, threads) for case in inputs for mode in ("score", "trace")
           for method in ("split", "split-all", "edlib", "wfa2", "biwfa")
           for threads in ((1, 2, 4, 8) if method == "split" else ((8,) if method == "split-all" else (1,)))]
fields = "block case mode method threads m n distance seconds peak_rss_kib proposed accepted nvpd_leaves wavefront_leaves retried".split()
scores = {}
# A measured call and its warm-up share the same wall-time budget.
# A censored configuration is recorded separately and removed from subsequent blocks.
censored = set()
with output.open("x") as stream:
    writer = csv.writer(stream)
    writer.writerow(fields)
    for block in range(1, 6):
        random.Random(8713 + block).shuffle(configs)
        for case, mode, method, threads in configs:
            key = case, mode, method, threads
            if key in censored:
                continue
            # taskset establishes affinity before libgomp discovers the placement.
            cpus = ",".join(str(c) for c in affinity[:threads])
            cmd = ["taskset", "-c", cpus, binary, inputs[case], method, mode, str(threads)]
            try:
                sample = subprocess.check_output(cmd, text=True, timeout=120).strip().split(",")
            except subprocess.TimeoutExpired:
                censored.add(key)
                dest = output.with_name("censored.csv")
                new = not dest.exists()
                with dest.open("a") as f:
                    failed = csv.writer(f)
                    if new:
                        failed.writerow(("block","case","mode","method","threads","process_limit_s"))
                    failed.writerow((block,*key,120))
                continue
            if len(sample) != 10 or scores.setdefault(case, sample[2]) != sample[2]:
                raise ValueError((case, sample))
            writer.writerow((block, case, mode, method, threads, *sample))
            stream.flush()
