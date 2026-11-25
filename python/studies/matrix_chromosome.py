"""Temps et efficacité d'une matrice distribuée de fenêtres chromosomiques."""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, INK, column_titles, estimate, label_panels, legend_above, use_style

ALGORITHMS = ("nvpd-omp", "nvpd-indexed", "edlib", "wfa2")
NAMES = ("NvPD", "NvPD indexé", "Edlib", "WFA2")


def read_trials(root):
    raw = pd.read_csv(root / "results/mpi/chromosome/839322/trials.csv")
    keys = ["records", "nodes", "algorithm", "schedule", "length_bound"]
    blocks = raw.groupby(keys).trial.agg(sorted)
    if len(raw) != 85 or len(blocks) != 17 or not all(b == [1, 2, 3, 4, 5] for b in blocks):
        raise ValueError("Chaque configuration exige cinq processus MPI distincts.")
    workload = raw.groupby("records")[["matrix_sha256", "equivalent_cells"]].nunique(dropna=False)
    if not workload.eq(1).all().all() or raw[["matrix_sha256", "equivalent_cells"]].isna().any().any():
        raise ValueError("La matrice et le nombre de cases doivent être identiques pour tous les solveurs.")
    if not ((raw.pairs == raw.records * (raw.records - 1) // 2) & (raw.cores == 192 * raw.nodes)).all():
        raise ValueError("Nombre de paires ou de cœurs incompatible avec le protocole.")
    return raw


def measures(raw):
    base = raw.loc[raw.schedule.eq("hilbert-weighted") & raw.length_bound.eq(0)]
    seconds = base.pivot(index=["algorithm", "trial"], columns=["records", "nodes"], values="seconds")
    area = base.groupby("records").equivalent_cells.first()
    strong = seconds[768, 1] / seconds[768, 2]
    weak = area[768] / (2 * area[544]) * seconds[544, 1] / seconds[768, 2]
    return pd.DataFrame({"seconds": seconds[768, 2], "speedup": strong, "strong": 50 * strong,
                         "weak": 100 * weak}).reset_index()


def make_figure(raw):
    use_style()
    data = measures(raw)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.65), sharey=True,
                             gridspec_kw={"width_ratios": (1.25, 1)})
    timing, efficiency = axes
    positions = np.arange(len(NAMES))[::-1]
    timing.set_yticks(positions, NAMES)
    for i, (key, pos) in enumerate(zip(ALGORITHMS, positions)):
        group = data.loc[data.algorithm.eq(key)]
        duration = sample(group.seconds)
        has_bound = raw.algorithm.eq(key) & raw.length_bound.eq(1)
        row = pos + .13 if has_bound.any() else pos
        timing.barh(row, duration["mean"], .22, color=COLORS[i], edgecolor=INK, linewidth=.6)
        estimate(timing, row, duration, INK, horizontal=True, marker="|")
        if has_bound.any():
            bounded = sample(raw.loc[has_bound, "seconds"])
            timing.barh(pos-.15, bounded["mean"], .2, facecolor="none", edgecolor=INK, linewidth=.65,
                        hatch="//", label="Avec test de sous-séquence" if i == 1 else None)
            estimate(timing, pos-.15, bounded, INK, horizontal=True, marker="|")
        fixed, scaled = sample(group.strong), sample(group.weak)
        efficiency.plot([fixed["mean"], scaled["mean"]], [pos+.10, pos-.10], color=INK, lw=.85)
        estimate(efficiency, pos+.10, fixed, COLORS[i], horizontal=True, marker="o",
                 label="Charge fixe" if i == 0 else None)
        estimate(efficiency, pos-.10, scaled, INK, horizontal=True, marker="s", hollow=True,
                 label="Charge par cœur fixée" if i == 0 else None)
    timing.set_xlim(left=0)
    timing.set_xticks([0,.4,.8,1.2,1.6])
    timing.set_xlabel("Temps sur 384 cœurs [s]")
    efficiency.set_xlabel("Efficacité de 192 à 384 cœurs [%]")
    efficiency.margins(x=.12)
    timing.set_ylim(-.55, len(NAMES)-.5)
    label_panels(axes)
    column_titles(axes, ("Durée de la matrice", "Efficacité parallèle"))
    legend_above(fig, axes, 3)
    return fig
