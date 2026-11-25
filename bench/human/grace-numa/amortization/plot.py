"""Durée des lots et amortissement d’une copie par puce Grace."""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "python"))
import ct_plot

DESCRIPTION = "Gain de temps de la copie locale sur deux puces Grace selon la taille du lot et durée du processus."


def direct_gains(summary, metric):
    rows = summary.loc[(summary.phase == "direct") & (summary.metric == metric)].sort_values("k")
    values = rows[["mean", "sd", "family_low95", "family_high95"]]
    if rows.empty or rows.k.duplicated().any() or not np.isfinite(values.to_numpy()).all():
        raise ValueError(f"Expected paired direct timings with simultaneous intervals for {metric}")
    return rows.set_index("k")


def near_crossing(gains, thresholds):
    """Show the dense sweep around the mean crossing and its simultaneous bounds."""
    ks = gains.index.to_numpy()
    runs = np.split(ks, np.flatnonzero(np.diff(ks) > 1) + 1)
    target = thresholds.loc["process_gain_s", "sustained_mean_k"]
    if pd.isna(target):
        target = ks[np.argmin(abs(gains["mean"]))]
    dense = [run for run in runs if len(run) >= 3] or runs
    run = min(dense, key=lambda part: np.min(abs(part - target)))
    landmarks = thresholds[["first_mean_k", "sustained_family_gain_k"]].stack()
    lower, upper = min(landmarks, default=target) - 5, max(landmarks, default=target) + 3
    zoom = run[(run >= lower) & (run <= upper)]
    return zoom if len(zoom) >= 4 else run


def draw_gain(axis, summary, thresholds):
    bounds = thresholds.loc[(thresholds.phase == "direct") & (thresholds.parity == "all")].set_index("metric")
    ks = near_crossing(direct_gains(summary, "process_gain_s"), bounds)
    axis.axhline(0, **ct_plot.GUIDE)
    for metric, color, marker, offset, label in (
        ("gain_s", ct_plot.COLORS[0], "o", -.13, "FASTA au score"),
        ("process_gain_s", ct_plot.COLORS[1], "s", .13, "Processus entier"),
    ):
        rows = direct_gains(summary, metric).loc[ks]
        center = rows["mean"].to_numpy() * 1000
        limits = rows[["family_low95", "family_high95"]].to_numpy().T * 1000
        axis.errorbar(ks + offset, center, yerr=rows.sd.to_numpy() * 1000,
                      fmt="none", ecolor=color, elinewidth=2.8, alpha=.75, zorder=2)
        axis.errorbar(ks + offset, center, yerr=np.stack([center - limits[0], limits[1] - center]),
                      fmt=marker, color=color, mec=ct_plot.INK, mew=.5, ecolor=color,
                      elinewidth=.85, capsize=2.3, ms=4, label=label, zorder=3)
    axis.set(xlabel="Nombre d’alignements [1]", ylabel="Gain de la copie locale [ms]",
             xlim=(ks.min() - .6, ks.max() + .6))
    axis.xaxis.set_major_locator(MaxNLocator(nbins=9, integer=True))
    axis.margins(y=.12)
    first, sustained = bounds.loc["process_gain_s", ["first_mean_k", "sustained_family_gain_k"]]
    note = []
    if pd.notna(first):
        note.append(f"Premier gain moyen à $K={first:.0f}$")
    if pd.notna(sustained):
        note.append(f"Borne positive dès $K={sustained:.0f}$")
    if note:
        axis.text(.025, .97, "Processus entier\n" + "\n".join(note),
                  transform=axis.transAxes, va="top", fontsize=8,
                  bbox=dict(facecolor=ct_plot.PAPER, edgecolor="none", pad=1.2))


def draw_totals(axis, summary):
    """Split process lifetime into FASTA input, replication, scoring and remaining time."""
    direct = summary.loc[summary.phase.eq("direct")].set_index(["k", "metric"])
    ks = direct.index.get_level_values("k")
    labels = []
    upper = 0.
    for group, k in enumerate((ks.min(), ks.max())):
        for j, mode, name in ((1, "local", "Local"), (0, "shared", "Partagé")):
            y = 1.3 * group + j * .48
            labels.append((y, f"{name} · {k}"))
            total = direct.loc[(k, f"{mode}_process_s")]
            read = direct.loc[(k, f"{mode}_read_s"), "mean"]
            copy = 0. if mode == "shared" else direct.loc[(k, "local_copy_s"), "mean"]
            align = direct.loc[(k, f"{mode}_align_s"), "mean"]
            residual = total["mean"] - read - copy - align
            left = 0.
            for width, color, label in ((read, ct_plot.COLORS[4], "Lecture FASTA"),
                                        (copy, ct_plot.COLORS[3], "Copie"),
                                        (align, ct_plot.COLORS[2], "Alignements"),
                                        (residual, ct_plot.COLORS[5], "Autres coûts")):
                if width < -1e-7:
                    raise ValueError("Measured components exceed process lifetime")
                if width > 0:
                    axis.barh(y, width, left=left, height=.36, color=color, edgecolor=ct_plot.INK,
                              linewidth=.45, label=label)
                left += width
            center, sd = total["mean"], total["sd"]
            half = total["high95"] - center
            axis.errorbar(center, y, xerr=sd, fmt="none", ecolor=ct_plot.INK, elinewidth=2)
            axis.errorbar(center, y, xerr=half, fmt="none", ecolor=ct_plot.INK, elinewidth=.7, capsize=2)
            endpoint = center + max(sd, half)
            upper = max(upper, endpoint)
            axis.annotate(f"{center:.2f}".replace(".", ","), (endpoint, y), xytext=(3, 0),
                          textcoords="offset points", va="center", fontsize=8,
                          bbox=dict(facecolor=ct_plot.PAPER, edgecolor="none", pad=.35))
    labels.sort()
    axis.set(yticks=[pos for pos, _ in labels], yticklabels=[name for _, name in labels],
             ylim=(labels[-1][0] + .4, -.4), xlim=(0, upper * 1.17), xlabel="Durée du processus [s]")
    axis.grid(False, axis="y")


def figure(summary, thresholds):
    ct_plot.use_style()
    fig, (gain, totals) = plt.subplots(1, 2, figsize=(8.8, 2.8), width_ratios=(1.4, 1))
    draw_gain(gain, summary, thresholds)
    draw_totals(totals, summary)
    ct_plot.label_panels((gain, totals))
    ct_plot.column_titles((gain, totals), ("Gain au voisinage du seuil", "Durée des lots courts et longs"))
    ct_plot.legend_above(fig, [gain, totals], 6)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    extension = args.output.suffix.lower()
    if extension not in {".png", ".pdf"}:
        parser.error("L’extension de sortie doit être .png ou .pdf")
    graph = figure(pd.read_csv(args.directory / "summary.csv"), pd.read_csv(args.directory / "threshold.csv"))
    with ct_plot.atomic_figure(args.output) as stream, plt.rc_context({"pdf.fonttype": 42}):
        if extension == ".pdf":
            graph.savefig(stream, format="pdf", metadata={"Title": DESCRIPTION, "CreationDate": None})
        else:
            graph.savefig(stream, format="png")
    plt.close(graph)
