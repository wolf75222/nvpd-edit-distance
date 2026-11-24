"""Sensibilité des ordonnanceurs MPI à l'ordre des fragments humains."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from ct_plot import use_style
from studies.measurement_matrix import measurement_matrix


def make_figure(root: Path):
    data = root / "results/mpi/schedules/839194"
    stats = pd.read_csv(data / "summary.csv")
    scaling = pd.read_csv(data / "scaling.csv")
    results = pd.concat([stats[stats.nodes.eq(2)], scaling]).set_index(["layout", "schedule", "metric"])
    policies = [("static", "Blocs"), ("cyclic", "Cyclique"), ("hilbert", "Hilbert"),
                ("hilbert-weighted", "Hilbert pondéré"), ("steal", "Vol")]
    layouts = ["locus", "length"]
    metrics = [("gcups", "Débit [GCUPS]", 1), ("max_mean_rank_cells", "Charge max./moy. [1]", 1),
               ("busy_fraction", "Occupation [%]", 100), ("efficiency", "Efficacité forte [%]", 100)]
    use_style()
    fig, axes = plt.subplots(1, 4, figsize=(9.0, 3.05))
    for j, (axis, (metric, title, factor)) in enumerate(zip(axes, metrics)):
        records = {(i, k): results.loc[(layout, policy, metric), ["mean", "sd", "high95"]] * factor
                   for i, (policy, _) in enumerate(policies) for k, layout in enumerate(layouts)}
        imbalance = metric == "max_mean_rank_cells"
        labels = [label if j == 0 else "" for _, label in policies]
        digits = ".2f" if imbalance else ".1f"
        measurement_matrix(axis, records, labels, ["Loci", "Longueurs"], reverse=imbalance,
                           precision=digits, spread=digits)
        axis.xaxis.tick_top()
        axis.set_xlabel(title, fontweight="bold", labelpad=6, fontsize=10)
        axis.xaxis.set_label_position("top")
    return fig
