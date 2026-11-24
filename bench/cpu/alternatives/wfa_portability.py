"""Check all exact scores and time paired stock/portable WFA2 in eight process blocks."""
import csv
import hashlib
import math
from pathlib import Path
import random
import statistics
import subprocess
import sys

build, inputs, out = map(Path, sys.argv[1:])
with (out / "machine-code.csv").open() as f:
    hashes = {(r["variant"], r["object"], r["part"]): r["sha256"] for r in csv.DictReader(f)}
if all(hashes["stock", name, "ELF"] == hashes["portable", name, "ELF"] for name in ("editdist", "warm")):
    print("Production and benchmark executables are byte-identical.")
    sys.exit(0)
with Path("results/alternatives/romeo/838762/cases.csv").open() as f:
    cases = [r for r in csv.DictReader(f) if r["case"].startswith(("20k-", "chr21-"))]
if len(cases) != 5:
    raise ValueError("expected five historical benchmark cases")
for c in cases:
    f = inputs / (c["case"] + ".fa")
    if hashlib.sha256(f.read_bytes()).hexdigest() != c["sha256"]:
        raise ValueError(f"input digest changed: {f}")
    check = subprocess.check_output([build/"stock/warm", f, "edlib", "1"], text=True)
    if int(check.split(",")[1]) != int(c["exact"]):
        raise ValueError("edlib reference changed")

order = random.Random(694621)
blocks = {}
with (out / "trials.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow("case variant block position call score seconds sha256".split())
    for block in range(1, 9):
        indices = list(range(len(cases)))
        order.shuffle(indices)
        for i in indices:
            case = cases[i]
            name, exact = case["case"], int(case["exact"])
            variants = ("stock", "portable") if (block+i)%2 else ("portable", "stock")
            for pos, variant in enumerate(variants, 1):
                output = subprocess.check_output([build/variant/"warm", inputs/(name+".fa"), "wfa2", "10"], text=True)
                runs = list(csv.reader(output.splitlines()))
                if len(runs) != 10:
                    raise ValueError(name)
                times = []
                for call, (seq, score, sec, _, _) in enumerate(runs, 1):
                    if int(seq) != call or int(score) != exact or float(sec) <= 0:
                        raise ValueError((case, call, score, sec))
                    times.append(float(sec))
                    w.writerow((name, variant, block, pos, call, score, sec, case["sha256"]))
                blocks[name, variant, block] = statistics.mean(times)
with (out / "summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow("case series units mean sd low95 high95".split())
    for c in cases:
        name = c["case"]
        stock = [blocks[name, "stock", b] for b in range(1,9)]
        portable = [blocks[name, "portable", b] for b in range(1,9)]
        ratios = [p/s for p,s in zip(portable,stock)]
        for series, values in (("stock_seconds", stock), ("portable_seconds", portable),
                               ("portable_over_stock", ratios)):
            mean, sd = statistics.mean(values), statistics.stdev(values)
            half = 2.364624251 * sd / math.sqrt(len(values))
            w.writerow((name, series, len(values), mean, sd, mean-half, mean+half))
