"""Compute paired strong and weak scaling on four complete CHM13 chromosome pairs."""
import argparse
from itertools import permutations
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "python"))
import ct_plot


def read(directory):
    directory = Path(directory)
    inputs = pd.read_csv(directory / "inputs.tsv", sep="\t", header=None,
                         names=["path", "m", "n", "score"])
    if [Path(p).stem for p in inputs.path] != ["chr1", "chr2", "chr3", "chr4"]:
        raise ValueError("The panel requires four distinct chromosome pairs in order")
    panel_bases = int((inputs.m + inputs.n).sum())
    rows = pd.read_csv(directory / "trials.csv", dtype={"cpu_ids": str})
    if (not np.isfinite(rows.seconds).all() or rows.seconds.le(0).any()
            or set(rows.block) != set(range(1, 7)) or not rows.workers.isin([1, 2, 4]).all()
            or not rows.phase.isin(["warmup", "measure"]).all()
            or rows.duplicated(["phase", "block", "position"]).any()
            or rows.duplicated(["phase", "block", "mode", "workers"]).any()):
        raise ValueError("Invalid or duplicate observation")
    warmup = rows.phase.eq("warmup")
    for block in rows.block.unique():
        block_rows = rows.loc[rows.block.eq(block)]
        plan = set(zip(block_rows.loc[~warmup, "mode"], block_rows.loc[~warmup, "workers"]))
        if (set(block_rows.loc[warmup, "workers"]) != {1, 2, 4}
                or not block_rows.loc[warmup, "mode"].eq("weak").all()
                or plan != {(mode, p) for mode in ("strong", "weak") for p in (1, 2, 4)}):
            raise ValueError(f"Incomplete configuration set in block {block}")
    passes = np.where(warmup, rows.workers, np.where(rows["mode"].eq("strong"), 16, 4 * rows.workers))
    if not rows.pairs.eq(4 * passes).all() or not rows.bases.eq(panel_bases * passes).all():
        raise ValueError("Logical workload disagrees with the input manifest")
    for row in rows.itertuples(index=False):
        if len(set(row.cpu_ids.split("|"))) != row.workers:
            raise ValueError("Thread placement disagrees with worker count")
    data = rows.loc[~warmup].copy()
    orders = {tuple(group.sort_values("position").workers.drop_duplicates())
              for _, group in data.groupby("block")}
    if orders != set(permutations([1, 2, 4])) or not rows.groupby("workers").cpu_ids.nunique().eq(1).all():
        raise ValueError("Unexpected block order or thread placement")
    return data


def analyze(data):
    paired = data.merge(data.loc[data.workers.eq(1), ["block", "mode", "seconds"]]
                        .rename(columns={"seconds": "baseline_s"}),
                        on=["block", "mode"], validate="many_to_one")
    ratio = paired.baseline_s / paired.seconds
    fixed = paired["mode"].eq("strong")
    paired["speedup"] = ratio.where(fixed)
    paired["efficiency"] = np.where(fixed, ratio / paired.workers, ratio)
    paired["throughput"] = paired.bases / paired.seconds / 1e9
    measures = paired.melt(id_vars=["mode", "workers"],
                           value_vars=["seconds", "speedup", "efficiency", "throughput"],
                           var_name="metric", value_name="value").dropna()
    summary = measures.groupby(["mode", "workers", "metric"], as_index=False).agg(
        n=("value", "count"), mean=("value", "mean"), sd=("value", "std"))
    half = t.ppf(.975, summary.n - 1) * summary.sd / np.sqrt(summary.n)
    return paired, summary.assign(low95=summary["mean"] - half, high95=summary["mean"] + half)


def decimal(value, precision):
    return f"{value:.{precision}f}".replace(".", ",")


def figure(summary):
    ct_plot.use_style()
    values = summary.set_index(["mode", "workers", "metric"])
    fig, ax = plt.subplots(figsize=(6.15, 2.65))
    for row, workers in enumerate((1, 2, 4)):
        for mode, offset, color, label in (("strong", -.19, ct_plot.COLORS[0], "Charge totale fixe"),
                                           ("weak", .19, ct_plot.COLORS[2], "Charge par cœur fixe")):
            y = row + offset
            seconds = values.loc[mode, workers, "seconds"]
            efficiency = values.loc[mode, workers, "efficiency"]["mean"]
            end = seconds["mean"]
            ax.barh(y, end, height=.30, color=color, edgecolor=ct_plot.INK, linewidth=.5,
                    label=label if row == 0 else None)
            for half, width, cap in ((seconds.high95 - end, .8, 3), (seconds.sd, 2.2, 0)):
                ax.errorbar(end, y, xerr=half, fmt="none", color=ct_plot.INK, lw=width, capsize=cap)
            ax.text(end / 2, y, f"$E$ = {decimal(100 * efficiency, 1)} %", ha="center", va="center",
                    fontsize=8)
            ax.annotate(decimal(end, 3), (seconds.high95, y), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=8, bbox=dict(facecolor=ct_plot.PAPER, edgecolor="none", pad=.3))
    ax.set(yticks=range(3), yticklabels=["1 cœur", "2 cœurs", "4 cœurs"],
           ylim=(2.48, -.48), xlim=(0, values.xs("seconds", level="metric").high95.max() * 1.12),
           xlabel="Temps du lot [s]")
    ax.xaxis.set_major_formatter(lambda x, _: decimal(x, 1))
    ax.grid(False, axis="y")
    ct_plot.legend_above(fig, [ax], 2)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    paired, summary = analyze(read(args.directory))
    paired.to_csv(args.directory / "paired.csv", index=False)
    summary.to_csv(args.directory / "summary.csv", index=False)
    if args.figure:
        graph = figure(summary)
        with ct_plot.atomic_figure(args.figure) as output:
            graph.savefig(output, format="png", dpi=220)
        plt.close(graph)


if __name__ == "__main__":
    main()
