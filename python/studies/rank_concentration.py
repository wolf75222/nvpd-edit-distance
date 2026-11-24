"""Répartition du travail dense entre les rangs d’une matrice de distances."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, GUIDE, INK, legend_above, use_style


def concentration(raw):
    """Part cumulée des cellules des k rangs les plus chargés de chaque lancement."""
    selected = raw.loc[raw.nodes.eq(2) & raw.layout.eq("length")]
    required = {"hilbert", "hilbert-weighted", "steal"}
    blocks = selected.groupby("schedule").block.agg(set).to_dict()
    if any(blocks.get(policy) != {1, 2, 3, 4, 5} for policy in required):
        raise ValueError("Les trois politiques doivent couvrir les cinq blocs de la campagne.")
    result = []
    for (policy, block), group in selected.groupby(["schedule", "block"]):
        if len(group) != 48 or group["rank"].nunique() != 48:
            raise ValueError("La répartition exige 48 rangs distincts par lancement.")
        if group.pairs.sum() != 8128 or group.cells.sum() != 181418301166:
            raise ValueError("Le triangle des 128 fragments est incomplet.")
        work = np.sort(group.cells.to_numpy())[::-1]
        share = np.r_[0., np.cumsum(work) / work.sum()]
        result.append(pd.DataFrame({"policy": policy, "block": block,
                                    "ranks": np.arange(len(share)), "share": share}))
    return pd.concat(result, ignore_index=True)


def make_figure(curves):
    use_style()
    fig, ax = plt.subplots(figsize=(6.1, 3.2))
    ax.plot([0, 48], [0, 100], **(GUIDE | {"zorder": 2.2, "lw": .6}), label="Partage égal")
    styles = (("hilbert", "Hilbert", COLORS[1], "o"), ("steal", "Vol", COLORS[0], "s"),
              ("hilbert-weighted", "Hilbert pondéré", COLORS[2], "D"))
    for policy, label, color, marker in styles:
        table = curves[curves.policy.eq(policy)].pivot(index="ranks", columns="block", values="share") * 100
        mean = table.mean(axis=1)
        if policy == "steal":
            intervals = [sample(row) for row in table.to_numpy()]
            ax.fill_between(table.index, [s["low95"] for s in intervals], [s["high95"] for s in intervals],
                            facecolor=color, alpha=.23, linewidth=0)
            ax.errorbar(12, mean.loc[12], yerr=intervals[12]["sd"], fmt="none",
                        ecolor=INK, capsize=3, linewidth=1.2)
        ax.plot(table.index, mean, color=color, label=label, linewidth=1.8)
        y = mean.loc[12]
        ax.plot(12, y, marker, color=color, markeredgecolor=INK, markersize=5)
        ax.annotate(f"{y:.1f} %".replace(".", ","), (12, y), xytext=(8, -13),
                    textcoords="offset points", fontsize=9)
    ax.set(xlabel="Effectif cumulé par charge décroissante [rangs]",
           ylabel="Part cumulée des cellules [%]", xticks=[0, 12, 24, 36, 48],
           yticks=[0, 25, 50, 75, 100], xlim=(0, 48), ylim=(0, 100))
    legend_above(fig, [ax], 4)
    return fig
