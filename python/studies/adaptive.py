"""Travail de la récurrence NvPD, temps et modèle mémoire empirique sur EPYC."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ct_plot import COLORS, GUIDE, INK, column_titles, estimate, label_panels, use_style
from studies.measurement_matrix import measurement_matrix


def make_figure(root: Path):
    """Comparer les sept paires et les étalons mesurés sur le même cœur EPYC."""
    data = root / "results/adaptive/839188"

    def load(name):
        with (data / name).open() as source:
            return list(csv.DictReader(source))

    def uncertainty(row, factor=1):
        return {k: float(row[k]) * factor for k in ("mean", "sd", "low95", "high95")}

    runs = {(r["case"], r["algorithm"], r["metric"]): r for r in load("summary.csv")}
    calib = {(r["kernel"], int(r["footprint_kib"])): r for r in load("calibration-summary.csv")}
    pmu = {(r["case"], r["algorithm"], r["metric"]): r for r in load("pmu-summary.csv")}
    work = {r["case"]: r for r in load("work.csv")}
    cases = ("syn-2", "syn-10", "independent", "chr21-25k", "chr21-100k", "chr21-250k", "chr21-1000k")
    labels = ("Synth. 2 %", "Synth. 10 %", "Indép.", "Chr21 25 k", "Chr21 100 k", "Chr21 250 k", "Chr21 1 M")
    methods = (("nvpd-adaptive", "Adaptatif"), ("nvpd-omp", "Dense"), ("banded-trim", "Bande"),
               ("edlib", "edlib"), ("wfa2", "WFA2"), ("biwfa", "BiWFA"))
    use_style()
    fig = plt.figure(figsize=(8.4, 5.85), layout="constrained")
    top, bottom = fig.add_gridspec(2, 1, height_ratios=(1.32, 1))
    time, ratio = top.subgridspec(1, 2, width_ratios=(2.1, 1)).subplots()
    bandwidth, roof = bottom.subgridspec(1, 2, width_ratios=(.9, 1.4)).subplots()

    times = {}
    for i, case in enumerate(cases):
        for j, (alg, _) in enumerate(methods):
            row = runs.get((case, alg, "seconds"))
            if row:
                times[i, j] = uncertainty(row, 1000)
    measurement_matrix(time, times, labels, [name for _, name in methods],
                       precision=".3g", spread=".1g", log=True)
    time.axhline(2.5, color=INK, linewidth=.8)
    column_titles((time, ratio), ("Temps par alignement [ms]", "Réduction du travail et des temps"))
    label_panels((time, ratio, bandwidth, roof))
    for i, case in enumerate(cases):
        ratio.plot(float(work[case]["cell_ratio"]), i-.2, marker="^", ls="none", color=COLORS[3],
                   markeredgecolor=INK, label="cellules denses / visitées" if i == 0 else None)
        for j, (alg, name) in enumerate((("nvpd-omp", "dense / adaptatif"),
                                          ("band-control", "scalaire / SIMD"))):
            row = runs.get((case, alg, "paired_time_over_adaptive"))
            if row:
                estimate(ratio, i+.19*j, uncertainty(row), COLORS[j], marker="os"[j],
                         horizontal=True, label=name if i == 0 else None)
    ratio.set(xscale="log", xlabel="rapport [1]", ylim=(6.5, -.5),
              yticks=range(7), yticklabels=[])
    ratio.axvline(1, **GUIDE)
    ratio.axhline(2.5, color=INK, linewidth=.8)
    fig.legend(*ratio.get_legend_handles_labels(), loc="outside upper center", ncol=3,
               fontsize=9, columnspacing=1.2, handlelength=1.3)

    sizes = (12, 384, 12288, 196608)
    for j, (kind, name) in enumerate((("read", "lecture"), ("triad", "triade"))):
        for i, size in enumerate(sizes):
            record = uncertainty(calib[kind, size])
            pos = i + (j-.5) * .37
            bandwidth.bar(pos, record["mean"], width=.34, color=COLORS[j], alpha=.57,
                          edgecolor=INK, linewidth=.6, label=name if i == 0 else None)
            estimate(bandwidth, pos, record, COLORS[j], marker="os"[j])
    bandwidth.set_xticks(range(4), ("12 Kio", "384 Kio", "12 Mio", "192 Mio"))
    bandwidth.set(xlabel="taille cumulée des tableaux", ylabel="débit utile [Go/s]", ylim=(0, None))
    bandwidth.legend(loc="upper right", fontsize=8, ncol=2, columnspacing=.8, handlelength=1)

    peak = uncertainty(calib["row", 12])
    read = uncertainty(calib["read", 196608])
    knee = peak["mean"] / read["mean"]
    x = np.geomspace(knee/2.5, 800, 300)
    roof.fill_between(x, np.minimum(peak["mean"]-peak["sd"], (read["mean"]-read["sd"])*x),
                      np.minimum(peak["mean"]+peak["sd"], (read["mean"]+read["sd"])*x),
                      color=INK, alpha=.11, linewidth=0)
    roof.plot(x, np.minimum(peak["mean"], read["mean"]*x), **GUIDE)
    roof.annotate(r"$P_{\rm L1}$", (knee*3, peak["mean"]), xytext=(0, 5),
                  textcoords="offset points", ha="center")
    for j, (case, name) in enumerate((("syn-10", "synth. 10 %"), ("independent", "indép."),
                                      ("chr21-100k", "chr21 100 k"))):
        ix = uncertainty(pmu[case, "nvpd-adaptive", "updates_per_estimated_dram_byte"])
        iy = uncertainty(pmu[case, "nvpd-adaptive", "gcell_s"])
        estimate(roof, iy["mean"], ix, COLORS[j], marker="os^"[j], horizontal=True)
        estimate(roof, ix["mean"], iy, COLORS[j], marker="os^"[j], label=name)
    roof.set(xscale="log", xlim=(x[0], x[-1]), ylim=(.8, 2.85),
             xlabel="intensité de lecture estimée [cellules/octet]", ylabel="débit [Gcell/s]")
    roof.legend(loc="lower right", ncol=2, fontsize=8, handlelength=1, columnspacing=.7)
    column_titles((bandwidth, roof), ("Étalonnage des tableaux", "Modèle de débit en lecture"))
    for axis in (ratio, roof):
        axis.xaxis.set_minor_locator(plt.NullLocator())
    return fig
