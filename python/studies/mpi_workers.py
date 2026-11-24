"""Effet apparié du moteur local sur la durée des lots de paires."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from ct_plot import COLORS, GUIDE, category_rows, column_titles, estimate, label_panels, legend_above, use_style


def make_figure(root: Path):
    data = pd.read_csv(root / "results/mpi/workers/839197/summary.csv")
    relative = data.loc[data.metric.eq("speedup_over_openmp")].set_index(["schedule", "nodes"])
    indexed = pd.read_csv(root / "results/mpi/chromosome/workers/839351/summary.csv")
    indexed = indexed.loc[indexed.metric.eq("openmp_over_kokkos")].set_index(["schedule", "nodes"])
    use_style()
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 5.25), height_ratios=[2.4, 1, 1.35])
    groups = [
        [("static", "Blocs"), ("cyclic", "Cyclique"), ("hilbert", "Hilbert"),
         ("hilbert-weighted", "Hilbert pondéré")],
        [("steal", "Vol")],
        [("hilbert-weighted", "Hilbert pondéré"), ("steal", "Vol")],
    ]
    for axis, group, results in zip(axes, groups, (relative, relative, indexed)):
        for i, (policy, _) in enumerate(group):
            for j, nodes in enumerate((1, 2)):
                row = results.loc[policy, nodes]
                stat = {key: (row[key] - 1) * 100 for key in ("mean", "high95")}
                stat["sd"] = 100 * row.sd
                estimate(axis, i + (j - .5) * .35, stat, COLORS[j], marker="os"[j],
                         label=f"{nodes * 192} cœurs", horizontal=True)
        axis.axvline(0, **GUIDE)
        category_rows(axis, [label for _, label in group], "")
    axes[-1].set_xlabel(r"Écart relatif $100\,(T_{\mathrm{OpenMP}}/T_{\mathrm{Kokkos}}-1)$ [%]")
    column_titles(axes, ["Affectation fixe · 8 128 paires", "Vol · 8 128 paires", "NvPD indexé · 294 528 paires"])
    label_panels(axes)
    legend_above(fig, axes, ncol=2)
    return fig
