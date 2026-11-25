"""Tracer les durées des alignements CHM13 et les écarts appariés entre placements."""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))
import ct_plot
from bench_stats import sample


MODES = {
    "single": "Un cœur, pages locales",
    "shared0": "Copie partagée sur A",
    "shared1": "Copie partagée sur B",
    "twins0": "Deux copies sur A",
    "twins1": "Deux copies sur B",
    "local": "Copies locales",
    "remote": "Copies distantes",
}
CONTRASTS = {
    "remote_minus_local": ("remote", "local", "Distant − local"),
    "twins0_minus_local": ("twins0", "local", "Copies sur A − locales"),
    "shared0_minus_twins0": ("shared0", "twins0", "Partagé A − copies sur A"),
    "shared0_minus_single": ("shared0", "single", "Partagé A − un cœur"),
}


def block_times(directory):
    trials = pd.read_csv(directory / "trials.csv")
    times = trials.pivot(index="block", columns="mode", values="seconds").sort_index()
    if (set(times.columns) != set(MODES) or len(times) < 2 or
            not np.isfinite(times.to_numpy()).all() or not (times > 0).all().all() or
            not trials.scores.eq(64).all() or trials.bases.nunique() != 1):
        raise ValueError("Every block must cover all placements at the same fixed load")
    return times


def delta_ms(times):
    return pd.DataFrame({key: 1000 * (times[left] - times[right])
                         for key, (left, right, _) in CONTRASTS.items()})


def figure(times):
    ct_plot.use_style()
    fig, (elapsed, delta) = plt.subplots(1, 2, figsize=(8.8, 3.25),
                                         width_ratios=(1.1, 1))
    ink, colors = ct_plot.INK, ct_plot.COLORS
    label_box = dict(facecolor=ct_plot.PAPER, edgecolor="none", pad=.5)
    upper = 0
    for y, mode in enumerate(MODES):
        stats = sample(times[mode])
        mean = stats["mean"]
        color = colors[4] if mode == "single" else colors[0]
        elapsed.barh(y, mean, height=.52, color=color, edgecolor=ink, linewidth=.55)
        elapsed.errorbar(mean, y, xerr=stats["sd"], fmt="none", ecolor=ink, elinewidth=2)
        elapsed.errorbar(mean, y, xerr=stats["high95"] - mean, fmt="none", ecolor=ink,
                         elinewidth=.8, capsize=2.5)
        bound = mean + max(stats["sd"], stats["high95"] - mean)
        upper = max(upper, bound)
        elapsed.annotate(f"{mean:.3f}".replace(".", ","), (bound, y), xytext=(3, 0),
                         textcoords="offset points", va="center", fontsize=8, bbox=label_box)
    elapsed.set(yticks=range(len(MODES)), yticklabels=list(MODES.values()),
                ylim=(len(MODES) - .43, -.57), xlim=(0, upper * 1.2), xlabel="Temps du lot [s]")
    elapsed.xaxis.set_major_formatter(lambda x, _: f"{x:.1f}".replace(".", ","))
    differences = delta_ms(times)
    delta.axvline(0, **ct_plot.GUIDE)
    for y, key in enumerate(CONTRASTS):
        values = differences[key]
        delta.hlines(y, 0, values.mean(), color=colors[2], linewidth=2.2, zorder=1.2)
        shifts = np.linspace(-.14, .14, len(values))
        delta.scatter(values, y + shifts, s=14, facecolors=ct_plot.PAPER,
                      edgecolors=colors[4], linewidths=.65, zorder=1.7,
                      label="Bloc" if y == 0 else None)
        stats = sample(values)
        ct_plot.estimate(delta, y, stats, colors[2], marker="D",
                         label="Différence moyenne" if y == 0 else None, horizontal=True)
        bound = max(values.max(), stats["high95"], stats["mean"] + stats["sd"])
        delta.annotate(f'{stats["mean"]:+.1f}'.replace(".", ","), (bound, y), xytext=(4, 0),
                       textcoords="offset points", va="center", fontsize=8, bbox=label_box)
    delta.set(yticks=range(len(CONTRASTS)), yticklabels=[v[2] for v in CONTRASTS.values()],
              ylim=(len(CONTRASTS) - .48, -.52), xlabel="Écart apparié [ms]")
    left, right = delta.dataLim.intervalx
    span = max(right - left, 1.)
    delta.set_xlim(left - .07 * span, right + .27 * span)
    for ax in (elapsed, delta):
        ax.grid(False, axis="y")
    ct_plot.label_panels((elapsed, delta))
    ct_plot.column_titles((elapsed, delta), ("Durée des 64 alignements", "Effet du placement"))
    ct_plot.legend_above(fig, [delta], 2)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    graph = figure(block_times(args.directory))
    with ct_plot.atomic_figure(args.output) as stream:
        graph.savefig(stream, format="png")
    plt.close(graph)
