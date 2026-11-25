"""Scores et accélérations pour les révisions chromosomiques de HG002."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

from bench_stats import sample
from ct_plot import COLORS, INK, column_titles, estimate, label_panels, legend_above, use_style


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


def read_flanks(root: Path):
    trials = pd.read_csv(root / "results/human/chromosome-dispatch/flanks/trials.csv")
    front = trials.loc[(trials.block > 0) & trials.method.isin(["sparse", "reuse"])]
    times = front.pivot(index=["chrom", "threads", "block"], columns="method", values="seconds")
    if (len(times) != 56 or times.isna().any().any() or (times <= 0).any().any()
            or not times.groupby(level=["chrom", "threads"]).size().eq(7).all()):
        raise ValueError("expected seven paired warm timings for each chromosome and core count")
    return (times["sparse"] / times["reuse"]).rename("flank_speedup")


def make_figure(trials, flanks):
    use_style()
    fig, (ax, after) = plt.subplots(1, 2, sharey=True, figsize=(8.75, 3.6),
                                    width_ratios=(1.15, 1))
    config = [("edlib", "Edlib"), ("split", "Fragments"), ("sparse", "Front borné")]
    chroms = ["chr1_MATERNAL", "chr1_PATERNAL", "chr21_MATERNAL", "chr21_PATERNAL"]
    labels = []
    limits = []
    for row, chrom in enumerate(chroms):
        case = trials[(trials.chrom == chrom) & (trials.threads == 1)]
        number, hap = chrom[3:].split("_")
        hap = "maternel" if hap == "MATERNAL" else "paternel"
        labels.append(f"Chr. {number} {hap}\n$d={int(case.score.iloc[0])}$")
        for k, (method, label) in enumerate(config):
            stat = sample(case.loc[case.method == method, "speedup_wfa"])
            limits.extend([stat["low95"], stat["high95"]])
            y = row + (k - 1) * .21
            ax.barh(y, abs(stat["mean"] - 1), left=min(1, stat["mean"]), height=.15,
                    color=COLORS[[4, 0, 2][k]], edgecolor=INK, linewidth=.5,
                    label=label if row == 0 else None)
            estimate(ax, y, stat, INK, marker="|", horizontal=True)
        for cores, offset, color, marker in [(1, -.13, COLORS[0], "o"), (8, .13, COLORS[3], "s")]:
            ratios = flanks.xs((chrom, cores), level=("chrom", "threads"))
            y = row + offset
            after.hlines(y, ratios.min(), ratios.max(), colors=INK, linewidth=.55)
            estimate(after, y, sample(ratios, independent=False), color, marker=marker, horizontal=True,
                     label=(f"{cores} cœur" if cores == 1 else f"{cores} cœurs") if row == 0 else None)
    ax.set_xscale("log")
    ax.set_xlim(min(x for x in limits if x > 0) * .8, max(limits) * 1.22)
    ticks = [x for x in (.01,.02,.05,.1,.2,.5,1,2,5,10) if ax.get_xlim()[0] <= x <= ax.get_xlim()[1]]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{x:g}".replace(".", ",") for x in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(range(len(chroms)), labels)
    ax.set_ylim(len(chroms) - .55, -.55)
    for axis in (ax, after):
        axis.grid(False, axis="y")
        for row in range(len(chroms)-1):
            axis.axhline(row + .5, color=INK, linestyle=":", linewidth=.5)
        axis.axvline(1, color=INK, linewidth=.85)
    ax.set_xlabel(r"$T_{\mathrm{WFA2},1} / T_{\mathrm{method},1}$ [1]")
    after.set_xlabel(r"$T_{\mathrm{rescan},p} / T_{\mathrm{reuse},p}$ [1]")
    after.margins(x=.09)
    label_panels([ax, after])
    column_titles([ax, after], ["Choix du solveur sur un cœur", "Réutilisation des extrémités"])
    legend_above(fig, [ax, after], ncol=5)
    return fig
