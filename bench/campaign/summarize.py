#!/usr/bin/env python3
"""Validate exact distances across runs and report per-pair median performance."""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

KEY = ("job", "host", "target", "binary_sha256", "input_sha256", "algorithm", "threads")
META = ("job", "host", "arch", "target", "binary_sha256", "input_sha256", "path", "length",
        "rate", "rep", "seed", "m", "n", "algorithm", "threads", "distance")
METRICS = ("trials", "median_seconds", "min_seconds", "max_seconds", "median_wall_seconds",
           "max_peak_rss_kib", "max_device_bytes", "equivalent_gcups", "dense_gcups",
           "sequential_reference_jobs", "speedup_vs_cpu_seq", "omp_efficiency_vs_cpu_seq", "omp_scaling", "omp_efficiency")
DENSE = {"seq", "diagonal", "nvpd-omp", "nvpd-cuda"}
OMP = {"diagonal", "nvpd-omp"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", nargs="+", type=Path)
    p.add_argument("--plan", type=Path, help="TSV with job and expected_trials; check coverage before summarizing")
    a = p.parse_args()
    job_counts = Counter()
    groups = defaultdict(list)
    distances = {}
    cpu_reference = defaultdict(list)
    skipped = 0
    seen = set()
    for path in a.files:
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                job_counts[row["job"]] += 1
                if row["status"] != "ok":
                    skipped += 1
                    continue
                seq = row["input_sha256"]
                reference = distances.setdefault(seq, (row["m"], row["n"], row["distance"]))
                if reference != (row["m"], row["n"], row["distance"]):
                    p.exit(2, f"Exact results disagree for {row['path']} in {path}\n")
                key = tuple(row[k] for k in KEY)
                trial = (*key, row["trial"])
                if trial in seen:
                    p.exit(2, f"Duplicate trial for {row['path']} in {path}\n")
                seen.add(trial)
                groups[key].append(row)
                if row["target"] == "cpu" and row["algorithm"] == "seq":
                    cpu_reference[seq].append(row)
    if a.plan:
        with a.plan.open(newline="") as stream:
            plan = {r["job"]: int(r["expected_trials"]) for r in csv.DictReader(stream, delimiter="\t")}
        for job in sorted(plan.keys() | job_counts.keys()):
            if job_counts[job] != plan.get(job):
                p.exit(2, f"Job {job}: observed {job_counts[job]} trials, expected {plan.get(job)}\n")
    if not groups:
        p.exit(2, "No successful measurements\n")

    def seconds(rows):
        return [float(x["seconds"]) for x in rows]

    out = csv.DictWriter(sys.stdout, fieldnames=META + METRICS, lineterminator="\n")
    out.writeheader()
    for key, rows in groups.items():
        row = rows[0]
        data = {k: row[k] for k in META}
        times = seconds(rows)
        t = median(times)
        equiv = int(row["m"]) * int(row["n"]) / t / 1e9
        data.update(trials=len(rows), median_seconds=t, min_seconds=min(times), max_seconds=max(times),
                    median_wall_seconds=median(float(x["wall_seconds"]) for x in rows),
                    max_peak_rss_kib=max(int(x["peak_rss_kib"]) for x in rows),
                    max_device_bytes=max(int(x["device_bytes"]) for x in rows), equivalent_gcups=equiv)
        if row["algorithm"] in DENSE:
            data["dense_gcups"] = equiv
        own_reference = groups.get((*key[:-2], "seq", "1")) if row["target"] == "cpu" else None
        reference_rows = own_reference or cpu_reference[row["input_sha256"]]
        if reference_rows:
            data["sequential_reference_jobs"] = "|".join(sorted({r["job"] for r in reference_rows}))
            speedup = median(seconds(reference_rows)) / t
            data["speedup_vs_cpu_seq"] = speedup
            if row["algorithm"] in OMP and row["target"] == "cpu":
                data["omp_efficiency_vs_cpu_seq"] = speedup / int(row["threads"])
        one = groups.get((*key[:-1], "1"))
        if one and row["algorithm"] in OMP:
            data["omp_scaling"] = median(seconds(one)) / t
            data["omp_efficiency"] = data["omp_scaling"] / int(row["threads"])
        out.writerow(data)
    if skipped:
        p.exit(1, f"{skipped} unsuccessful trials appear in the raw CSV files.\n")


if __name__ == "__main__":
    main()
