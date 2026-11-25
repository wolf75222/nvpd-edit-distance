#!/usr/bin/env python3
"""Summarize timing, residency and observed dispatch over independent 100-knt pairs."""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

PARENT = Path(__file__).resolve().parents[1]
sys.path.append(str(PARENT))
from focus_analysis import ct_plot, summarize

METHODS = ("nvpd-cuda", "nvpd-cuda-persistent", "nvpd-hybrid-score", "wfa2", "edlib")
LABELS = {"nvpd-cuda": "NvPD CUDA", "nvpd-cuda-persistent": "NvPD persistant",
          "nvpd-hybrid-score": "Choix CPU/GPU", "edlib": "Edlib"}
GROUP_LABELS = {"100k-2pct": "2 %", "control-5pct": "5 %", "original-five": "10 %",
                "held-out": "10 % suppl."}


def read(directory):
    raw = pd.read_csv(directory / "trials.csv")
    manifest = pd.read_csv(directory / "inputs.tsv", sep="\t")
    keys = ["group", "rep"]
    if manifest.duplicated(keys).any() or manifest.sha256.duplicated().any():
        raise ValueError("Duplicate input in manifest")
    groups = set(manifest.group)
    required = {"100k-2pct", "control-5pct", "original-five"}
    if not required <= groups <= GROUP_LABELS.keys() or not manifest.groupby("group").size().eq(5).all():
        raise ValueError("Expected five pairs at each declared divergence")
    data = raw.merge(manifest[keys + ["sha256", "distance"]], on=keys, how="left",
                     suffixes=("", "_oracle"), validate="many_to_one")
    if (not data.distance.eq(data.distance_oracle).all() or data.sha256.nunique() != len(manifest)
            or set(data.method) != set(METHODS) or set(data.block) != set(range(1, 6)) or set(data.order) != set(range(1, 6))
            or data.duplicated(keys + ["block", "method"]).any() or data.duplicated(keys + ["block", "order"]).any()
            or not data.groupby(keys + ["block"]).method.nunique().eq(5).all()
            or not data.groupby(keys + ["method"]).order.nunique().eq(5).all()
            or not data.groupby(keys + ["method"]).size().eq(5).all()):
        raise ValueError("Incomplete or unbalanced five-block design")
    hybrid = data.method.eq("nvpd-hybrid-score")
    if (not data.groupby(keys + ["method"]).backend.nunique().eq(1).all() or
            not data.loc[hybrid, "backend"].isin(["wfa2", "nvpd-cuda", "nvpd-cuda-persistent"]).all()
            or not data.loc[~hybrid, "backend"].eq(data.loc[~hybrid, "method"]).all()):
        raise ValueError("Unexpected hybrid dispatch")
    measurements = data[["seconds", "warmup_seconds", "after_warmup_pss_kib", "after_warmup_rss_kib"]]
    if (not np.isfinite(measurements).all().all() or measurements.le(0).any().any()
            or not data.after_warmup_rss_kib.ge(data.after_warmup_pss_kib).all()
            or data.device_bytes.lt(0).any() or not data.loc[data.backend.isin(["edlib", "wfa2"]), "device_bytes"].eq(0).all()):
        raise ValueError("Invalid time or memory observation")
    return data


def analyze(data):
    index = ["group", "rep", "block"]
    times = data.pivot(index=index, columns="method", values="seconds")
    comparisons = []
    for method, reference in ([(m, "wfa2") for m in LABELS]
                              + [("nvpd-cuda-persistent", "nvpd-cuda"), ("nvpd-hybrid-score", "nvpd-cuda")]):
        comparisons.append((times[method] / times[reference]).rename("ratio").reset_index()
                           .assign(method=method, reference=reference))
    selected = data.query("method == 'nvpd-hybrid-score'").set_index(index).backend
    underlying = times.stack().reindex(
        pd.MultiIndex.from_frame(selected.reset_index().rename(columns={"backend": "method"}))
    ).set_axis(selected.index)
    comparisons.append((times["nvpd-hybrid-score"] / underlying).rename("ratio").reset_index()
                       .assign(method="nvpd-hybrid-score", reference="selected"))
    all_ratios = pd.concat(comparisons, ignore_index=True)
    ratio_keys = ["group", "method", "reference"]
    pairs = summarize(all_ratios.groupby(ratio_keys + ["rep"]).ratio)
    wins = pairs.assign(win=pairs["mean"].lt(1)).groupby(ratio_keys).win.sum()
    groups = summarize(pairs.groupby(ratio_keys)["mean"]).merge(wins.rename("wins"), on=ratio_keys)
    resource_keys = ["group", "method", "metric"]
    values = data.melt(id_vars=["group", "method", "rep"], var_name="metric", value_name="value",
                       value_vars=["seconds", "after_warmup_rss_kib", "after_warmup_pss_kib", "device_bytes"])
    pair_values = values.groupby(resource_keys + ["rep"]).value.mean()
    resources = summarize(pair_values.groupby(level=resource_keys))
    routes = data.query("method == 'nvpd-hybrid-score'")[["group", "rep", "sha256", "backend"]].drop_duplicates()
    return pairs, groups, resources, routes


def number(x, digits):
    return f"{x:.{digits}f}".replace(".", ",")


def figure(summary):
    ct_plot.use_style()
    ratios = summary.query("reference == 'wfa2'").set_index(["group", "method"])
    cohorts = [g for g in GROUP_LABELS if g in ratios.index.get_level_values("group")]
    fig, ax = plt.subplots(figsize=(2.55 + 1.42 * len(cohorts), 2.8))
    colors = (ct_plot.COLORS[2], ct_plot.COLORS[4], ct_plot.COLORS[1])
    names = ("Plus rapide que WFA2", "IC contenant 1", "Plus lent que WFA2")
    seen = set()
    for row, (method, name) in enumerate(LABELS.items()):
        ax.text(-.07, row + .5, name, ha="right", va="center")
        for col, group in enumerate(cohorts):
            r = ratios.loc[(group, method)]
            category = 0 if r.high95 < 1 else 2 if r.low95 > 1 else 1
            ax.add_patch(Rectangle((col + .015, row + .025), .97, .95, color=colors[category],
                                   alpha=.37, linewidth=0, label=None if category in seen else names[category]))
            seen.add(category)
            digits = 3 if r["mean"] < 2 else 2
            ax.text(col+.5, row+.32, f"{number(r['mean'], digits)} ± {number(r.sd, digits)}",
                    ha="center", va="center", fontweight="bold")
            ax.text(col+.5, row+.69, f"[{number(r.low95, digits)} ; {number(r.high95, digits)}]",
                    ha="center", va="center", fontsize=8)
    for col, group in enumerate(cohorts):
        title = GROUP_LABELS[group] if group == "held-out" else GROUP_LABELS[group] + " d’éditions"
        ax.text(col+.5, -.16, title, ha="center", va="bottom", fontweight="bold")
    ax.text(-.07, -.16, r"$T/T_{\mathrm{WFA2}}$ [1]", ha="right", va="bottom", fontweight="bold")
    ax.set(xlim=(-.72, len(cohorts)+.01), ylim=(len(LABELS)+.01, -.55))
    ax.set_axis_off()
    ct_plot.legend_above(fig, [ax], ncol=len(seen))
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timing", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    frames = analyze(read(args.timing))
    args.output.mkdir(parents=True, exist_ok=True)
    for name, frame in zip(("pairs", "summary", "resources", "routes"), frames):
        frame.to_csv(args.output / f"{name}.csv", index=False)
    if args.figure:
        graph = figure(frames[1])
        graph.savefig(args.figure)
        plt.close(graph)


if __name__ == "__main__":
    main()
