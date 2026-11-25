#!/usr/bin/env python3
"""Confirm score routing on independent pairs. Run the workers or summarize raw observations."""
import argparse
from collections import defaultdict
import csv
from decimal import Decimal
import json
from math import isfinite, sqrt
import os
from pathlib import Path
import random
from statistics import mean, stdev
import subprocess
from tempfile import TemporaryDirectory

from measure import measure, pin_timing_process, sha256
from synthetic import generate

RATES = {"0.02": "2pct", "0.05": "5pct", "0.1": "10pct"}
REPS = range(12, 17)
OLD, NEW = "hybrid-original", "hybrid-refined"
HYBRID = "nvpd-hybrid-score"
METHODS = (OLD, NEW, "nvpd-cuda", "nvpd-cuda-persistent", "wfa2", "edlib")
ROUTES = {OLD: {"nvpd-cuda", "wfa2"}, NEW: {"nvpd-cuda", "nvpd-cuda-persistent", "wfa2"}}
METRICS = ("seconds", "after_warmup_rss_kib", "after_warmup_pss_kib", "device_bytes")


def write_rows(path, rows, delimiter=","):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], delimiter=delimiter, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_rows(path, delimiter=","):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def pair(cli, algorithm, fasta):
    text = subprocess.check_output([cli, "pair", "--algorithm", algorithm, "--threads", "1",
                                    "--input", fasta], text=True)
    return json.loads(text)


def prepare(directory, builds):
    generate(directory, [100000], list(map(Decimal, RATES)), max(REPS), 2026)
    items = []
    for row in read_rows(directory / "manifest.tsv", "\t"):
        if int(row["rep"]) not in REPS:
            continue
        fasta = directory / row["path"]
        distance = int(pair(builds[NEW] / "editdist", "edlib", fasta)["distance"])
        item = dict(group=RATES[row["rate"]], **row, sha256=sha256(fasta), distance=distance)
        for method, routes in ROUTES.items():
            cli = pair(builds[method] / "editdist", HYBRID, fasta)
            if int(cli["distance"]) != distance or cli["backend"] not in routes:
                raise ValueError(f"Invalid {method} CLI answer on {fasta}: {cli}")
            item[method] = cli["backend"]
        items.append(item)
    return items


def run(reference, candidate, output):
    output.mkdir(parents=True, exist_ok=False)
    pin_timing_process(output)
    builds = {OLD: reference.resolve(), NEW: candidate.resolve()}
    write_rows(output / "binaries.tsv", [
        dict(job=os.environ.get("SLURM_JOB_ID", "local"), build=name, program=exe, sha256=sha256(path / exe))
        for name, path in builds.items() for exe in ("repeat_kernel", "editdist")], "\t")
    rng = random.Random(410231)
    with TemporaryDirectory(prefix="refine-pairs-") as temp:
        directory = Path(temp)
        items = prepare(directory, builds)
        write_rows(output / "inputs.tsv", items, "\t")
        with (output / "trials.csv").open("x", newline="") as fh:
            writer = None
            for block in range(6):
                order = list(enumerate(items))
                rng.shuffle(order)
                for idx, item in order:
                    # Six Williams orders balance both position and first-order carryover per pair.
                    for position, offset in enumerate((0, 1, 5, 2, 4, 3), 1):
                        method = METHODS[(idx + block + offset) % 6]
                        binary = builds.get(method, builds[NEW]) / "repeat_kernel"
                        record = measure(binary, HYBRID if method in ROUTES else method, directory, item)
                        if method in ROUTES and record["backend"] != item[method]:
                            raise ValueError(f"CLI and worker routes differ for {item['path']}, {method}")
                        row = dict(group=item["group"], rep=item["rep"], block=block+1, order=position,
                                   method=method, distance=item["distance"], **record)
                        if writer is None:
                            writer = csv.DictWriter(fh, fieldnames=row, lineterminator="\n")
                            writer.writeheader()
                        writer.writerow(row)
                    fh.flush()
                print(f"Completed block {block+1}/6", flush=True)


def stats(values):
    # Two-sided Student intervals: df=5 across blocks or df=4 across independent pairs.
    crit = {5: 2.7764451051977987, 6: 2.570581835636314}[len(values)]
    mu, sd = mean(values), stdev(values)
    half = crit * sd / sqrt(len(values))
    return dict(count=len(values), mean=mu, sd=sd, low95=mu-half, high95=mu+half)


def analyze(timing, output):
    manifest = read_rows(timing / "inputs.tsv", "\t")
    items = {(m["group"], int(m["rep"])): m for m in manifest}
    expected = {(group, rep) for group in RATES.values() for rep in REPS}
    if (set(items) != expected or len(manifest) != len(expected)
            or len({r["sha256"] for r in manifest}) != len(expected)):
        raise ValueError("Missing or duplicated independent pairs")
    data = {}
    for row in read_rows(timing / "trials.csv"):
        group, rep, block, method = row["group"], int(row["rep"]), int(row["block"]), row["method"]
        key = (group, rep, block, method)
        item = items[group, rep]
        backend = row["backend"]
        if (key in data or method not in METHODS or not 1 <= block <= 6 or row["distance"] != item["distance"]
                or backend != (item[method] if method in ROUTES else method)
                or not all(isfinite(float(row[k])) and float(row[k]) > 0 for k in ("seconds", "warmup_seconds"))
                or not 0 < int(row["after_warmup_pss_kib"]) <= int(row["after_warmup_rss_kib"])
                or min(int(row["device_bytes"]), int(row["warmup_device_bytes"])) < 0
                or (backend in ("wfa2", "edlib") and int(row["device_bytes"]) != 0)):
            raise ValueError(f"Invalid trial {key}")
        data[key] = row
    pair_results, groups = [], defaultdict(list)
    for group, rep in sorted(items):
        item = items[group, rep]
        for block in range(1, 7):
            if {int(data[group, rep, block, m]["order"]) for m in METHODS} != set(range(1, 7)):
                raise ValueError("Duplicate position within a block")
        for method in METHODS:
            rows = [data[group, rep, block, method] for block in range(1, 7)]
            if {int(r["order"]) for r in rows} != set(range(1, 7)):
                raise ValueError("Unbalanced sequence of methods")
            series = {(metric, ""): [float(r[metric]) for r in rows] for metric in METRICS}
            references = [m for m in (METHODS if method == NEW else ("wfa2",)) if m != method]
            if method in ROUTES:
                references.append("selected")
            for reference in references:
                ref = item[method] if reference == "selected" else reference
                series["ratio", reference] = [float(r["seconds"]) / float(data[group, rep, block, ref]["seconds"])
                                              for block, r in enumerate(rows, 1)]
            for (metric, ref), values in series.items():
                result = stats(values)
                groups[group, method, metric, ref].append(result["mean"])
                pair_results.append(dict(group=group, rep=rep, sha256=item["sha256"], method=method,
                                         backend=rows[0]["backend"], metric=metric, reference=ref, **result))
    summary = [dict(group=g, method=m, metric=k, reference=r, **stats(vals),
                    wins=sum(v < 1 for v in vals) if k == "ratio" else "")
               for (g, m, k, r), vals in groups.items()]
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "pairs.csv", pair_results)
    write_rows(output / "summary.csv", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    run_args = sub.add_parser("run")
    for field in ("reference", "candidate", "output"):
        run_args.add_argument(field, type=Path)
    analyze_args = sub.add_parser("analyze")
    for field in ("timing", "output"):
        analyze_args.add_argument(field, type=Path)
    args = parser.parse_args()
    if args.action == "run":
        run(args.reference, args.candidate, args.output)
    else:
        analyze(args.timing, args.output)
