"""GH200 integer roof and end-to-end throughput on two chromosome regions."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, INK, column_titles, estimate, join, label_panels, legend_above

CALIBRATION_JOB = "839177"
LOCAL_JOB = "839191"


def make_figure(root: Path):
    data = root / "results" / "roofline" / "gpu"
    calibration = pd.read_csv(data / CALIBRATION_JOB / "calibration.csv")
    trials = pd.read_csv(data / LOCAL_JOB / "trials.csv")
    figure, axes = plt.subplots(1, 3, figsize=(9.6, 3.05), layout="constrained")
    left = axes[0]
    copy = calibration[calibration.method == "stream-copy"]
    reg = calibration[calibration.method == "register-updates"]
    bandwidth = (copy.requested_bytes / copy.seconds / 1e9).mean()
    ceiling = (reg.updates / reg.seconds / 1e9).mean()
    intensity = np.geomspace(1 / 22, 24, 200)
    left.plot(
        intensity,
        np.minimum(bandwidth * intensity, ceiling),
        color=INK,
        linestyle="--",
        linewidth=0.9,
        label="Enveloppe empirique",
    )
    for steps, rows in calibration[calibration.method == "stream-updates"].groupby("steps"):
        stats = sample(rows.updates / rows.seconds / 1e9)
        estimate(
            left,
            steps / 16,
            stats,
            COLORS[4],
            marker="D",
            label="Tableaux de 512 Mio" if steps == 1 else None,
        )
    left.set(
        xscale="log",
        yscale="log",
        xlabel="Intensité [mises à jour/octet demandé]",
        ylabel="Débit [Gmises à jour/s]",
        xlim=(1 / 22, 24),
        ylim=(110, 8300),
    )
    left.set_xticks([1 / 16, 1 / 4, 1, 4, 16], ["1/16", "1/4", "1", "4", "16"])
    left.set_yticks([200, 500, 1000, 2000, 5000], ["200", "500", "1 000", "2 000", "5 000"])

    styles = [
        ("native", "Une ligne", COLORS[0], "o"),
        ("fused", "Deux lignes", COLORS[1], "s"),
        ("local", "Halo local", COLORS[2], "^"),
    ]
    regions = ["collinear", "insertion"]
    lengths = [16384, 65536, 131072, 524288]
    for axis, region in zip(axes[1:], regions):
        for method, name, color, marker in styles:
            series = trials[(trials.method == method) & trials.case.str.startswith(region + "-")]
            points = []
            for n in lengths:
                rows = series[series.m == n]
                stats = sample(rows.m * rows.n / rows.seconds / 1e9)
                estimate(axis, n, stats, color, marker=marker, label=name if n == lengths[0] else None)
                points.append(stats["mean"])
            join(axis, lengths, points, color)
        axis.set(xscale="log", xlabel="Longueur GRCh38 [knt]", ylabel="Débit logique [GCUPS]")
        axis.set_xticks(lengths, ["16,4", "65,5", "131", "524"])
        axis.set_ylim(bottom=0)
    ymax = max(ax.get_ylim()[1] for ax in axes[1:])
    for ax in axes[1:]:
        ax.set_ylim(0, ymax)
    legend_above(figure, axes, ncol=5)
    column_titles(axes, ["Mise à jour entière", "Région colinéaire", "Région avec insertion"])
    label_panels(axes)
    return figure
