#!/usr/bin/env python3
"""Select by geometric mean of paired speedups on two 32k chromosome windows."""
import csv
import math
import statistics
import sys
from collections import defaultdict

with open(sys.argv[1], newline="") as f:
    measurements = list(csv.DictReader(f))
print("threads,selected,iterations")
for threads in ("1", "4"):
    cases = defaultdict(dict)
    near = defaultdict(list)
    for r in measurements:
        if r["length"] == "32768" and r["threads"] == threads:
            cases[r["cohort"], r["block"]][r["variant"]] = float(r["seconds"])
            if r["cohort"] == "near":
                near[r["variant"]].append(float(r["seconds"]))
    score = {v: statistics.mean(math.log(case["baseline"] / case[v]) for case in cases.values()) for v in near}
    selected = max(score, key=score.get)
    repeats = math.ceil(22 / min(statistics.median(near[v]) for v in ("baseline", selected)))
    print(f"{threads},{selected},{repeats}")
