#!/usr/bin/env python3
"""Compute blocked uncertainty and break-even workloads from independent fresh-process confirmation."""

import argparse
from collections import defaultdict
import csv
import functools
import math
from pathlib import Path
import random
from statistics import mean, stdev

from measure_e2e import MAX_K, MODES, PREFIX_BLOCKS, SEED, TIMES, balanced, confirm_plan


@functools.lru_cache(maxsize=None)
def student_quantile(df, prob):
    """Positive Student quantile, using the finite trigonometric recurrence for integer df."""
    def integral(theta):
        power = df - 1
        even = power % 2 == 0
        value = theta if even else math.sin(theta)
        for m in range(2 if even else 3, power + 1, 2):
            value = ((m - 1) * value + math.sin(theta) * math.cos(theta)**(m - 1)) / m
        return value
    norm = integral(math.pi / 2)
    lo, hi = 0., math.pi / 2
    for _ in range(60):
        mid = (lo + hi) / 2
        if .5 + .5 * integral(mid) / norm < prob:
            lo = mid
        else:
            hi = mid
    return math.sqrt(df) * math.tan((lo + hi) / 2)


def estimates(sample, family_size):
    n, m, sd = len(sample), mean(sample), stdev(sample)
    se = sd / math.sqrt(n)
    half = student_quantile(n - 1, .975) * se
    result = dict(n=n, mean=m, sd=sd, low95=m-half, high95=m+half,
                  family_low95="", family_high95="")
    if family_size:
        h = student_quantile(n - 1, 1 - .025 / family_size) * se
        result.update(family_low95=m-h, family_high95=m+h)
    return result


def csv_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify(data, direct_plan):
    prefix = [r for r in data if r["phase"] == "prefix"]
    cross, direct_keys, blocks, _ = confirm_plan(prefix)
    plan_k = [int(r["k"]) for r in direct_plan]
    if (plan_k != direct_keys or any(int(r["direct_blocks"]) != blocks for r in direct_plan)
            or any(r["pilot_sustained_k"] != (str(cross) if cross else "") for r in direct_plan)):
        raise ValueError("Direct workload selection differs from the calibration rule")
    keys = {"prefix": list(range(2, MAX_K + 1, 2)), "direct": direct_keys}
    block_counts = {"prefix": PREFIX_BLOCKS, "direct": blocks}
    indexed = {}
    for r in data:
        phase, mode = r["phase"], r["mode"]
        k, b, rotation, position = (int(r[s]) for s in ("k", "block", "rotation", "position"))
        group = phase, k, b, mode
        if (phase not in keys or mode not in MODES or k not in keys[phase] or b not in range(1, block_counts[phase]+1)
                or int(r["scores"]) != k or group in indexed):
            raise ValueError(f"Invalid or repeated trial: {group}")
        expected_rot, reverse = balanced(block_counts[phase], k if phase == "direct" else 0)[b-1]
        within = 1
        if phase == "direct":
            order = list(direct_keys)
            random.Random(SEED + 1000 + b).shuffle(order)
            within = order.index(k) + 1
        if rotation != expected_rot or (MODES[::-1] if reverse else MODES).index(mode) + 1 != position or int(r["order"]) != within:
            raise ValueError(f"Unmatched workload ordering/rotation in {group}")
        t = {field: float(r[field]) for field in TIMES}
        if phase == "direct":
            t["process_s"] = float(r["process_s"])
            if t["process_s"] < t["total_s"]:
                raise ValueError(f"Invalid process timer in {group}")
        elif r["process_s"]:
            raise ValueError(f"Prefix process time cannot be assigned to a partial workload: {group}")
        if (not all(math.isfinite(v) and v >= 0 for v in t.values()) or
                min(t["read_s"], t["align_s"]) <= 0 or (mode == "shared") != (t["copy_s"] == 0) or
                abs(t["read_s"] + t["copy_s"] - t["setup_s"]) > 1e-5 or
                abs(t["setup_s"] + t["align_s"] - t["total_s"]) > 1e-5):
            raise ValueError(f"Invalid timing decomposition in {group}")
        indexed[group] = t
    if len(indexed) != 2 * sum(len(ks) * block_counts[p] for p, ks in keys.items()):
        raise ValueError("Missing an expected timing")
    for mode in MODES:
        for b in range(1, PREFIX_BLOCKS + 1):
            partial = [indexed["prefix", k, b, mode] for k in keys["prefix"]]
            if (len({t["setup_s"] for t in partial}) != 1 or
                    any(a["total_s"] >= z["total_s"] for a, z in zip(partial, partial[1:]))):
                raise ValueError("Prefix observations did not share one start with advancing checkpoints")
    return indexed, keys, block_counts


def thresholds(results, phase, metric, parity, ks):
    ks = [k for k in ks if parity == "all" or k % 2 == (1 if parity == "odd" else 0)]
    if not ks:
        return None
    def first(cond):
        return next((k for k in ks if cond(results[k])), "")
    def sustained(key):
        return next((k for i, k in enumerate(ks) if all(results[x][key] > 0 for x in ks[i:])), "")
    return dict(phase=phase, metric=metric, parity=parity, tested_k="|".join(map(str, ks)),
                first_mean_k=first(lambda r: r["mean"] > 0), sustained_mean_k=sustained("mean"),
                first_possible_family_gain_k=first(lambda r: r["family_high95"] > 0),
                first_family_gain_k=first(lambda r: r["family_low95"] > 0),
                last_family_loss_k=next((k for k in reversed(ks) if results[k]["family_high95"] < 0), ""),
                sustained_family_gain_k=sustained("family_low95"))


def analyze(folder):
    trials, keys, block_counts = verify(csv_rows(folder/"timings.csv"), csv_rows(folder/"confirm-design.csv"))
    pairs, summary, bounds = [], [], []
    for phase, ks in keys.items():
        gains = defaultdict(dict)
        for k in ks:
            series = defaultdict(list)
            for b in range(1, block_counts[phase] + 1):
                a, z = (trials[phase, k, b, mode] for mode in MODES)
                metrics = {"gain_s": a["total_s"] - z["total_s"], "setup_cost_s": z["setup_s"] - a["setup_s"],
                           "align_gain_s": a["align_s"] - z["align_s"], "local_copy_s": z["copy_s"],
                           "read_difference_s": a["read_s"] - z["read_s"]}
                for mode, data in (("shared", a), ("local", z)):
                    for field in ("read_s", "setup_s", "align_s", "total_s"):
                        metrics[mode+"_"+field] = data[field]
                if phase == "direct":
                    metrics.update(shared_process_s=a["process_s"], local_process_s=z["process_s"],
                                   process_gain_s=a["process_s"] - z["process_s"])
                pairs.append(dict(phase=phase, k=k, block=b, **metrics))
                for name, value in metrics.items():
                    series[name].append(value)
            for name, values in series.items():
                contrast = name in ("gain_s", "process_gain_s")
                row = dict(phase=phase, k=k, metric=name, **estimates(values, len(ks) if contrast else 0))
                summary.append(row)
                if contrast:
                    gains[name][k] = row
        for name, results in gains.items():
            for parity in ("all", "even", "odd"):
                row = thresholds(results, phase, name, parity, ks)
                if row:
                    bounds.append(row)
    # Prefix checkpoints share a process and have no process-lifetime contrast.
    for row in pairs:
        for name in ("shared_process_s", "local_process_s", "process_gain_s"):
            row.setdefault(name, "")
    write_csv(folder/"paired.csv", pairs)
    write_csv(folder/"summary.csv", summary)
    write_csv(folder/"threshold.csv", bounds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    analyze(parser.parse_args().directory)
