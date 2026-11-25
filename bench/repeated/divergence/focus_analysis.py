#!/usr/bin/env python3
"""Summarize paired NvPD comparisons across independent sequence pairs."""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullLocator
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))
import ct_plot

METHODS = ("original", "candidate", "edlib", "wfa2")
RATIOS = (("original", "wfa2"), ("candidate", "wfa2"), ("candidate", "original"),
          ("original", "edlib"), ("candidate", "edlib"))
LABELS = {"original-five": "100 knt, 10 % · initiales", "held-out": "100 knt, 10 % · supplémentaires",
          "control-5pct": "100 knt, 5 %", "100k-2pct": "100 knt, 2 %", "1m-2pct": "1 Mnt, 2 %"}
KEYS = ["group", "rep", "block"]


def summarize(series):
    grouped = series.agg(["count", "mean", "std"]).rename(columns={"std": "sd"})
    half = t.ppf(.975, grouped["count"] - 1) * grouped.sd / np.sqrt(grouped["count"])
    return grouped.assign(low95=grouped["mean"] - half, high95=grouped["mean"] + half).reset_index()


def read_trials(directories):
    frames = []
    seen_groups = set()
    for directory in directories:
        trials = pd.read_csv(directory / "trials.csv")
        manifest = pd.read_csv(directory / "inputs.tsv", sep="\t")
        groups_here = set(manifest.group)
        if seen_groups & groups_here or not groups_here <= LABELS.keys():
            raise ValueError(f"Duplicate or unknown group in {directory}")
        seen_groups.update(groups_here)
        if (manifest.duplicated(["group", "rep"]).any() or manifest.sha256.duplicated().any()
                or not trials.groupby(["group", "rep"]).size().index.equals(
                    manifest.set_index(["group", "rep"]).sort_index().index)):
            raise ValueError(f"Invalid paired-input coverage in {directory}")
        trials = trials.merge(manifest[["group", "rep", "distance", "sha256"]], on=["group", "rep"],
                              suffixes=("", "_oracle"), validate="many_to_one")
        if not trials.distance.eq(trials.distance_oracle).all():
            raise ValueError(f"Score disagrees with input manifest in {directory}")
        frames.append(trials)
    data = pd.concat(frames, ignore_index=True)
    if data.duplicated(KEYS + ["method"]).any() or data.duplicated(KEYS + ["order"]).any():
        raise ValueError("Duplicate observation or method position in a block")
    if (set(data.method) != set(METHODS) or not data.groupby(KEYS).size().eq(4).all()
            or not data.order.isin([1, 2, 3, 4]).all()
            or not data.groupby(["group", "rep"]).sha256.nunique().eq(1).all()
            or not data.groupby("group").sha256.nunique().eq(5).all()):
        raise ValueError("Each group needs five distinct pairs and four methods per block")
    position_counts = data.groupby(["group", "rep", "method", "order"]).size().unstack("order")
    if (position_counts.isna().any().any() or not position_counts.nunique(axis=1).eq(1).all()
            or data.groupby(["group", "rep"]).block.nunique().groupby("group").nunique().ne(1).any()):
        raise ValueError("Unbalanced block design")
    times = data[["seconds", "warmup_seconds"]]
    if not np.isfinite(times).all().all() or (times <= 0).any().any():
        raise ValueError("Invalid alignment time")
    if data[["rss_kib", "device_bytes"]].isna().any().any() or (data[["rss_kib", "device_bytes"]] < 0).any().any():
        raise ValueError("Invalid host or device memory footprint")
    return data


def analyze(data):
    times = data.pivot(index=KEYS, columns="method", values="seconds")
    ratios = pd.concat(
        [(times[method] / times[reference]).rename("ratio").reset_index()
         .assign(method=method, reference=reference) for method, reference in RATIOS], ignore_index=True)
    ratio_keys = ["group", "method", "reference"]
    pairs = summarize(ratios.groupby(ratio_keys + ["rep"]).ratio)
    groups = pairs.groupby(ratio_keys)
    summary = summarize(groups["mean"]).merge(
        pairs.assign(win=pairs["mean"].lt(1)).groupby(ratio_keys).win.sum().rename("wins").reset_index(),
        on=ratio_keys)
    # Each pair has equal weight, including when pooling elapsed time and host/device peaks.
    # On Grace getrusage can report a zero peak even for a resident process.
    timing_resources = data.rename(columns={"rss_kib": "ru_maxrss_kib"}).copy()
    timing_resources["ru_maxrss_kib"] = timing_resources.ru_maxrss_kib.replace(0, np.nan)
    measurements = timing_resources.melt(id_vars=["group", "rep", "method"],
                             value_vars=["seconds", "ru_maxrss_kib", "device_bytes"],
                             var_name="metric", value_name="value")
    group_keys = ["group", "method", "metric"]
    pair_means = measurements.groupby(group_keys + ["rep"]).value.mean()
    resources = summarize(pair_means.groupby(level=group_keys))
    return pairs, summary, resources


def sampled_memory(path, data):
    samples = pd.read_csv(path)
    if (samples.duplicated(KEYS + ["method"]).any() or samples.sampled_smaps_rss_kib.le(0).any()
            or samples.samples.lt(5).any() or not samples.groupby(KEYS).method.nunique().eq(4).all()):
        raise ValueError("Incomplete paired RSS measurements")
    identities = data[["group", "rep", "distance", "sha256"]].drop_duplicates()
    check = samples.merge(identities, on=["group", "rep", "distance", "sha256"], how="left",
                          indicator=True, validate="many_to_one")
    if not check["_merge"].eq("both").all() or not samples.groupby(["group", "method"]).rep.nunique().eq(5).all():
        raise ValueError("Timing and RSS experiments refer to different pairs")
    pairs = samples.groupby(["group", "method", "rep"]).sampled_smaps_rss_kib.mean()
    return summarize(pairs.groupby(level=["group", "method"])).assign(metric="sampled_smaps_rss_kib")


def cohort_results(path):
    raw = pd.read_csv(path)
    regimes = {"shared", "concurrent"}
    if (set(raw["mode"]) != regimes or raw.duplicated(["mode", "block"]).any()
            or not raw.groupby("block").size().eq(2).all() or not raw.groupby("block").rep_order.nunique().eq(1).all()
            or raw.worker_sha256.nunique() != 1 or not raw.workers.eq(raw["mode"].map({"shared": 1, "concurrent": 5})).all()):
        raise ValueError("Invalid cohort execution modes, workloads or block matching")
    metrics = ["seconds", "sum_pss_kib", "sum_rss_kib", "sum_worker_peak_device_bytes"]
    values = raw.melt(id_vars=["block", "mode"], value_vars=metrics, var_name="metric")
    if not np.isfinite(values.value).all() or values.value.le(0).any():
        raise ValueError("Missing cohort measurement")
    result = summarize(values.groupby(["mode", "metric"]).value)
    paired = values.pivot(index=["metric", "block"], columns="mode", values="value")
    ratios = (paired.concurrent / paired.shared).rename("ratio")
    return pd.concat([result, summarize(ratios.groupby(level="metric")).assign(mode="concurrent/shared")],
                     ignore_index=True)



def figure(summary):
    ct_plot.use_style()
    methods = (("original", "NvPD initial", "s", ct_plot.COLORS[0]),
               ("candidate", "NvPD optimisé", "o", ct_plot.COLORS[2]))
    wfa = summary.loc[summary.reference.eq("wfa2")].set_index(["group", "method"])
    groups = [key for key in LABELS if key in wfa.index.get_level_values("group")]
    fig = plt.figure(figsize=(7.1, 1.25 + .48 * len(groups)))
    grid = fig.add_gridspec(1, 2, width_ratios=[5.9, 1])
    ax = fig.add_subplot(grid[0, 0])
    counts = fig.add_subplot(grid[0, 1], sharey=ax)
    counts.set_axis_off()
    counts.set_xlim(0, 1)
    counts.text(.5, 1.02, "Paires plus rapides", transform=counts.transAxes,
                ha="center", va="bottom", fontsize=8)
    for i, group in enumerate(groups):
        initial = wfa.loc[(group, "original"), "mean"]
        candidate = wfa.loc[(group, "candidate"), "mean"]
        ct_plot.join(ax, (i - .14, i + .14), (initial, candidate), ct_plot.COLORS[4], horizontal=True)
        digits = 3 if candidate < 2 else 2
        ax.annotate(f"{candidate:.{digits}f}".replace(".", ","), (candidate, i + .14),
                    xytext=(-6, 0), textcoords="offset points", ha="right", va="center", fontsize=8,
                    bbox={"facecolor": ct_plot.PAPER, "edgecolor": "none", "pad": .25})
        for (name, label, marker, color), delta in zip(methods, (-.14, .14)):
            record = wfa.loc[group, name]
            ct_plot.estimate(ax, i + delta, record, color, marker=marker,
                             label=label if i == 0 else None, horizontal=True)
            counts.plot(.27, i + delta, marker=marker, ms=3.4, linestyle="none", color=color,
                        markeredgecolor=ct_plot.INK, markeredgewidth=.45)
            counts.text(.41, i + delta, f"{int(record.wins)}/{int(record['count'])}",
                        va="center", fontsize=8)
    lower, upper = min(1.0, wfa.low95.min()), max(1.0, wfa.high95.max())
    if lower > 0 and upper / lower > 4:
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1, ls="--", lw=.8, color=ct_plot.INK)
    ct_plot.category_rows(ax, [LABELS[g] for g in groups], "Temps NvPD / temps WFA2 [1]")
    ct_plot.legend_above(fig, [ax], ncol=2)
    ax.xaxis.set_major_formatter(lambda x, _: f"{x:g}".replace(".", ","))
    return fig


def memory_figure(resources, cohorts):
    """GPU workspace peaks by input length and proportional host RSS by process layout."""
    ct_plot.use_style()
    fig, (gpu, host) = plt.subplots(1, 2, figsize=(8.1, 2.65), gridspec_kw={"width_ratios": [1.02, 1]})
    device = resources.loc[resources.metric.eq("device_bytes")].set_index(["group", "method"])
    fields = ["mean", "sd", "low95", "high95"]
    methods = (("original", "NvPD initial", "s", ct_plot.COLORS[0]),
               ("candidate", "NvPD optimisé", "o", ct_plot.COLORS[2]))
    lower, upper = float("inf"), 0
    for row, group in enumerate(("100k-2pct", "1m-2pct")):
        estimates = []
        for method, label, marker, color in methods:
            values = device.loc[(group, method), fields] / 2**20
            lower = min(lower, values.low95, values["mean"] - values.sd)
            upper = max(upper, values.high95, values["mean"] + values.sd)
            estimates.append((values, color, marker, label if row == 0 else None))
        ct_plot.comparison(gpu, row, estimates, digits=2)
    ct_plot.category_rows(gpu, ["100 knt", "1 Mnt"], "Pic des tampons GPU [Mio]")
    if lower > 0 and upper / lower > 4:
        gpu.set_xscale("log")
        gpu.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        gpu.xaxis.set_minor_locator(NullLocator())
        gpu.set_xlim(lower / 1.75, upper * 1.95)
    else:
        gpu.set_xlim(0, upper * 1.25)
    gpu.xaxis.set_major_formatter(lambda x, _: f"{x:g}".replace(".", ","))

    observed = cohorts.set_index(["mode", "metric"])
    tick_labels = []
    right = 0
    for row, (mode, label) in enumerate((("shared", "Un processus"), ("concurrent", "Cinq processus"))):
        pss = observed.loc[(mode, "sum_pss_kib"), fields] / 1024
        wall = observed.loc[(mode, "seconds"), fields]
        right = max(right, pss.high95)
        tick_labels.append(label + "\n" + f"{wall['mean']:.4f} ± {wall.sd:.4f} s".replace(".", ","))
        host.barh(row, pss["mean"], height=.27, color=ct_plot.COLORS[4], alpha=.45,
                  edgecolor=ct_plot.INK, linewidth=.6)
        ct_plot.estimate(host, row, pss, ct_plot.INK, marker="|", horizontal=True)
        ct_plot.value_label(host, pss.high95, row, pss["mean"], digits=1, horizontal=True)
    ct_plot.category_rows(host, tick_labels, "Somme des PSS hôte [Mio]")
    host.set_xlim(0, right * 1.19)
    ct_plot.column_titles((gpu, host), ("Alignement à 2 % d’éditions", "Lot de cinq paires à 10 %"))
    ct_plot.label_panels((gpu, host))
    ct_plot.legend_above(fig, (gpu,), ncol=2)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+", help="disjoint groups collected with focus.py")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--memory-figure", type=Path, help="workspace and cohort PSS figure")
    parser.add_argument("--smaps", type=Path, help="independent RSS observations from focus_memory.py")
    parser.add_argument("--cohort", type=Path, help="cohort memory and batch timings from cohort_memory.py")
    args = parser.parse_args()
    data = read_trials(args.inputs)
    pairs, summary, resources = analyze(data)
    if args.smaps:
        resources = pd.concat([resources, sampled_memory(args.smaps, data)], ignore_index=True)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, frame in [("pairs", pairs), ("summary", summary), ("resources", resources)]:
        frame.to_csv(args.output / f"{name}.csv", index=False)
    if args.memory_figure and not args.cohort:
        parser.error("--memory-figure requires --cohort")
    if args.cohort:
        cohort = cohort_results(args.cohort)
        cohort.to_csv(args.output / "cohort-summary.csv", index=False)
        if args.memory_figure:
            memory = memory_figure(resources, cohort)
            memory.savefig(args.memory_figure)
            plt.close(memory)
    if args.figure:
        plot = figure(summary)
        plot.savefig(args.figure)
        plt.close(plot)


if __name__ == "__main__":
    main()
