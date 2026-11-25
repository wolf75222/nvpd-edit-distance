#!/usr/bin/env python3
"""Estimate paired timing contrasts and audit the placement of both chromosome panels."""

import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import re
from statistics import mean, stdev

from measure import CPUS, MODES, design, input_lengths, panel_owners

T95 = 2.1603686564610127  # Two-sided Student quantile, df = 13.
BLOCKS = range(1, 15)
CONTRASTS = (("remote", "local"), ("shared0", "single"), ("shared1", "single"),
             ("single", "local"), ("single", "remote"), ("shared0", "local"), ("shared1", "local"),
             ("shared0", "twins0"), ("twins0", "local"), ("shared1", "twins1"), ("twins1", "local"),
             ("twins1", "twins0"))


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(sample):
    m, sd = mean(sample), stdev(sample)
    h = T95 * sd / math.sqrt(len(sample))
    return dict(count=len(sample), mean=m, sd=sd, low95=m - h, high95=m + h)


def verify_trials(directory, lengths):
    workload = 16 * sum(lengths.values())
    orders = design()
    values = {}
    for row in read_csv(directory / "trials.csv"):
        b, p = int(row["block"]), int(row["position"])
        mode = row["mode"]
        workers = 1 if mode == "single" else 2
        if (b not in BLOCKS or p not in range(1, 8) or mode != orders[b-1][p-1] or (b, mode) in values
                or int(row["workers"]) != workers or int(row["scores"]) != 64 or int(row["bases"]) != workload
                or row["cpu_ids"] != ("0" if workers == 1 else "0|72")):
            raise ValueError(f"Wrong trial, workload or core binding: {row}")
        seconds = float(row["seconds"])
        if not all(math.isfinite(x) and x > 0 for x in (seconds, float(row["warmup_s"]))):
            raise ValueError(f"Invalid duration: {row}")
        values[b, mode] = seconds
    if len(values) != len(BLOCKS) * len(MODES):
        raise ValueError("Missing a paired timing")
    return values


def verify_residency(directory, lengths):
    """Check all 98 × 2 snapshots against the VMA extent, expected node and panel owner."""
    topology = (directory / "topology.txt").read_text()
    page_size = int(re.search(r"^base_page_bytes=(\d+)$", topology, re.M)[1])
    nodes = {int(cpu): int(node) for cpu, node in re.findall(r"^cpu=(\d+) .*? node=(\d+) ", topology, re.M)}
    if set(nodes) != set(CPUS) or len(set(nodes.values())) != 2:
        raise ValueError(f"Incorrect NUMA topology: {nodes}")
    observations = {}
    for row in read_csv(directory / "residency.csv"):
        b, mode, phase = int(row["block"]), row["mode"], row["phase"]
        c, k, s = (int(row[f]) for f in ("copy", "chrom", "side"))
        owner, expected, node = (int(row[f]) for f in ("owner_cpu", "expected_node", "observed_node"))
        length, mapped, amount = (int(row[f]) for f in ("bytes", "mapped_bytes", "resident_bytes"))
        key = (b, mode, phase, c, k, s)
        if (b not in BLOCKS or mode not in MODES or phase not in ("before", "after")
                or c not in range(len(panel_owners(mode))) or owner != panel_owners(mode)[c] or expected != nodes[owner]
                or lengths.get((k, s)) != length or mapped != (length + page_size - 1) // page_size * page_size
                or amount <= 0 or amount % page_size or node < 0):
            raise ValueError(f"Invalid residency observation: {row}")
        counts = observations.setdefault(key, {})
        if node in counts:
            raise ValueError(f"Duplicate residency observation for node {node}: {key}")
        counts[node] = amount
    panels = defaultdict(list)
    for b in BLOCKS:
        for mode in MODES:
            for phase in ("before", "after"):
                for c, owner in enumerate(panel_owners(mode)):
                    owned = total = 0
                    for (k, s), length in lengths.items():
                        counts = observations.get((b, mode, phase, c, k, s), {})
                        mapped = (length + page_size - 1) // page_size * page_size
                        if sum(counts.values()) != mapped:
                            raise ValueError(f"Incomplete residency: {b}/{mode}/{phase}/{c}/{k}/{s}")
                        owned += counts.get(nodes[owner], 0)
                        total += mapped
                    panels[mode, phase, c, nodes[owner]].append((owned / total, total))
    rows = []
    for (mode, phase, c, node), sample in panels.items():
        ratios, sizes = zip(*sample)
        rows.append(dict(mode=mode, phase=phase, copy=c, expected_node=node, blocks=len(sample),
                         resident_mib=mean(sizes) / 2**20, min_expected_pct=100*min(ratios),
                         mean_expected_pct=100*mean(ratios), max_expected_pct=100*max(ratios)))
    write_csv(directory / "residency-summary.csv", rows)
    under = [r for r in rows if r["min_expected_pct"] <= 99]
    if under:
        raise ValueError(f"Expected >99% of each panel resident on its assigned node: {under}")


def fieller(numer, denom):
    """Ratio of paired means and Fieller 95% interval, if a positive excess is resolved."""
    n, a, b = len(numer), mean(numer), mean(denom)
    if summarize(denom)["low95"] <= 0:
        return dict(ratio="", low95="", high95="")
    cov = sum((x - a) * (y - b) for x, y in zip(numer, denom)) / (n - 1)
    t = T95 * T95 / n
    qa, qb, qc = b*b - t*stdev(denom)**2, -2*(a*b-t*cov), a*a - t*stdev(numer)**2
    disc = math.sqrt(max(0., qb*qb - 4*qa*qc))
    return dict(ratio=a/b, low95=(-qb-disc)/(2*qa), high95=(-qb+disc)/(2*qa))


def analyze(directory, manifest):
    lengths = input_lengths(manifest)
    values = verify_trials(directory, lengths)
    verify_residency(directory, lengths)
    pairs = []
    for b in BLOCKS:
        row = dict(block=b)
        for mode in MODES:
            row[f"time_{mode}_s"] = values[b, mode]
        for a, base in CONTRASTS:
            x, y = values[b, a], values[b, base]
            row[f"{a}_minus_{base}_s"] = x - y
            row[f"{a}_over_{base}"] = x / y
        pairs.append(row)
    write_csv(directory / "paired.csv", pairs)
    columns = {k: [row[k] for row in pairs] for k in pairs[0] if k != "block"}
    write_csv(directory / "summary.csv", [dict(metric=k, **summarize(v)) for k, v in columns.items()])
    # The denominator is the slowdown observed in the current matched experiment.
    denom = columns["shared0_minus_single_s"]
    shares = []
    for key in ("shared0_minus_twins0_s", "twins0_minus_local_s", "shared0_minus_local_s"):
        v = columns[key]
        shares.append(dict(numerator=key, denominator="shared0_minus_single_s", numerator_mean_s=mean(v),
                           denominator_mean_s=mean(denom), **fieller(v, denom)))
    write_csv(directory / "attribution.csv", shares)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    manifest = args.manifest or args.directory / "inputs.tsv"
    if not manifest.is_file() and args.manifest is None:
        manifest = Path(__file__).resolve().parents[3] / "results/human/grace-scaling/inputs.tsv"
    analyze(args.directory, manifest)
