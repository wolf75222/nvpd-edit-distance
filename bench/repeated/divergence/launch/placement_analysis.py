#!/usr/bin/env python3
"""Summarize the six paired placement blocks on the same Grace node."""
import csv
import math
from pathlib import Path
import statistics
import sys

MODES = ("master", "owner")
BLOCKS = range(1, 7)
T95 = 2.570581835636305  # Student t quantile at 0.975 with five degrees of freedom.


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def summarize_residency(directory):
    rows = read_csv(directory / "residency.csv")
    if not rows and (directory / "residency-unavailable.txt").exists():
        return
    nodes = sorted({int(row["numa_node"]) for row in rows})
    scopes = sorted({row["scope"] for row in rows})
    if "sequence_vmas" not in scopes:
        raise ValueError("Missing residency observation for the sequence buffers")
    amounts = {}
    for row in rows:
        key = (int(row["block"]), row["mode"], row["scope"], int(row["numa_node"]))
        amounts[key] = amounts.get(key, 0) + int(row["resident_pages"]) * int(row["page_kib"])
    with (directory / "residency-summary.csv").open("w", newline="") as stream:
        result = csv.writer(stream)
        result.writerow(("scope", "mode", "numa_node", "mean_MiB", "sd_MiB", "mean_percent", "sd_percent"))
        for scope in scopes:
            for mode in MODES:
                per_block = [[amounts.get((block, mode, scope, node), 0) / 1024 for node in nodes] for block in BLOCKS]
                if any(sum(row) <= 0 for row in per_block):
                    raise ValueError(f"Incomplete residency sample for {scope}/{mode}")
                for i, node in enumerate(nodes):
                    resident = [row[i] for row in per_block]
                    percent = [100 * row[i] / sum(row) for row in per_block]
                    result.writerow((scope, mode, node, statistics.mean(resident), statistics.stdev(resident),
                                     statistics.mean(percent), statistics.stdev(percent)))


def summarize(directory):
    runs = read_csv(directory / "trials.csv")
    keys = {(int(x["block"]), x["mode"]) for x in runs}
    if len(runs) != 12 or keys != {(b, mode) for b in BLOCKS for mode in MODES}:
        raise ValueError("Incomplete or duplicate placement blocks")
    first = [row["mode"] for row in runs if row["position"] == "1"]
    if len(first) != 6 or first.count("master") != 3 or any(
            {x["position"] for x in runs if int(x["block"]) == b} != {"1", "2"} for b in BLOCKS):
        raise ValueError("Unbalanced crossover order")
    placements = {row["cpu_ids"] for row in runs}
    if len(placements) != 1 or len(set(placements.pop().split("|"))) != 4:
        raise ValueError("Core placement changed")
    manifest = directory / "inputs.tsv"
    if not manifest.exists():  # Archived runs kept the manifest with the scaling trials.
        manifest = directory.parent / "inputs.tsv"
    with manifest.open() as stream:
        bases = 16 * sum(int(row[1]) + int(row[2]) for row in csv.reader(stream, delimiter="\t"))
    values = {}
    for row in runs:
        if int(row["pairs"]) != 64 or int(row["bases"]) != bases:
            raise ValueError("Unexpected logical workload")
        read, score = float(row["read_seconds"]), float(row["seconds"])
        if not all(math.isfinite(x) and x > 0 for x in (read, score, float(row["warmup_seconds"]))):
            raise ValueError("Invalid duration")
        values[int(row["block"]), row["mode"]] = {"read": read, "score": score, "read_plus_score": read + score}
    with (directory / "paired.csv").open("w", newline="") as paired, (directory / "summary.csv").open("w", newline="") as summary:
        details, aggregates = csv.writer(paired), csv.writer(summary)
        details.writerow(("block", "phase", "master_s", "owner_s", "master_over_owner", "master_minus_owner_s"))
        aggregates.writerow(("phase", "metric", "blocks", "mean", "sd", "low95", "high95"))
        for phase in ("read", "score", "read_plus_score"):
            columns = {"master_s": [], "owner_s": [], "master_over_owner": [], "master_minus_owner_s": []}
            for block in BLOCKS:
                master, owner = (values[block, mode][phase] for mode in MODES)
                observation = master, owner, master / owner, master - owner
                details.writerow((block, phase, *observation))
                for samples, value in zip(columns.values(), observation):
                    samples.append(value)
            for metric, samples in columns.items():
                mean, sd = statistics.mean(samples), statistics.stdev(samples)
                half = T95 * sd / math.sqrt(len(samples))
                aggregates.writerow((phase, metric, len(samples), mean, sd, mean - half, mean + half))


if __name__ == "__main__":
    directory = Path(sys.argv[1])
    summarize(directory)
    summarize_residency(directory)
