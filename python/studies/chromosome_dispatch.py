"""Scores et accélérations pour les révisions chromosomiques de HG002."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

from bench_stats import sample
from ct_plot import COLORS, INK, estimate, legend_above, use_style


def read_trials(root: Path):
    trials = pd.read_csv(root / "results/human/chromosome-dispatch/trials.csv")
    groups = ["chrom", "method", "threads"]
    if trials.duplicated(groups + ["block"]).any() or not (trials.groupby(groups).size() == 5).all():
        raise ValueError("expected five complete timing blocks per configuration")
    if (len(trials) != 220 or trials.chrom.nunique() != 4 or
            not (trials.groupby("chrom").score.nunique() == 1).all()):
        raise ValueError("incomplete chromosome experiment or inconsistent score")
    keys = ["chrom", "block"]
    baseline = trials.loc[trials.method == "wfa2", keys + ["seconds"]]
    trials = trials.merge(baseline.rename(columns={"seconds": "wfa_s"}), on=keys, validate="many_to_one")
    trials["speedup_wfa"] = trials.wfa_s / trials.seconds
    single = trials.loc[trials.threads == 1, keys + ["method", "seconds"]]
    trials = trials.merge(single.rename(columns={"seconds": "single_s"}), on=keys + ["method"],
                          validate="many_to_one")
    trials["speedup"] = trials.single_s / trials.seconds
    trials["efficiency"] = trials.speedup / trials.threads
    trials["rss_mib"] = trials.rss_kib / 1024
    return trials


def make_figure(trials):
    use_style()
    fig, ax = plt.subplots(figsize=(8.6, 3.55))
    config = [("edlib", 1, "edlib · 1"), ("biwfa", 1, "BiWFA · 1"),
              ("split", 8, "Fragments · 8"), ("sparse", 1, "Front · 1"),
              ("sparse", 8, "Front · 8")]
    chroms = ["chr1_MATERNAL", "chr1_PATERNAL", "chr21_MATERNAL", "chr21_PATERNAL"]
    labels = []
    limits = []
    palette = [COLORS[i] for i in (4, 1, 0, 3, 2)]
    for row, chrom in enumerate(chroms):
        case = trials[trials.chrom == chrom]
        number, hap = chrom[3:].split("_")
        hap = "maternel" if hap == "MATERNAL" else "paternel"
        labels.append(f"Chr. {number} {hap}\n$d={int(case.score.iloc[0])}$")
        for i, (method, threads, label) in enumerate(config):
            ratios = case[(case.method == method) & (case.threads == threads)].speedup_wfa
            stat = sample(ratios)
            limits.extend([stat["low95"], stat["mean"], stat["high95"]])
            y = row + (i - 2) * .15
            ax.barh(y, abs(stat["mean"] - 1), left=min(1, stat["mean"]), height=.115,
                    color=palette[i], edgecolor=INK, linewidth=.45,
                    label=label if not row else None, zorder=2)
            estimate(ax, y, stat, INK, marker="|", horizontal=True)
    ax.set_xscale("log")
    positive = [x for x in limits if np.isfinite(x) and x > 0]
    ax.set_xlim(min(positive) * .85, max(positive) * 1.16)
    ticks = [x for x in (.01,.02,.05,.1,.2,.5,1,2,5,10,20,50,100) if ax.get_xlim()[0] <= x <= ax.get_xlim()[1]]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{x:g}".replace(".", ",") for x in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(range(len(chroms)), labels)
    ax.set_ylim(len(chroms)-.55, -.55)
    ax.grid(False, axis="y")
    for row in range(len(chroms)-1):
        ax.axhline(row+.5, color=INK, linestyle=":", linewidth=.5)
    ax.axvline(1, color=INK, linewidth=1)
    ax.set_xlabel(r"Rapport des temps $T_{\mathrm{WFA2},1}/T_{a,p}$ [1]")
    legend_above(fig, [ax], len(config))
    return fig
