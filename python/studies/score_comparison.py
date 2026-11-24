"""Rapports appariés des temps NvPD, Edlib et WFA2 sur un même nœud GH200."""
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

from ct_plot import COLORS, INK, use_style


def paired_ratios(raw):
    data = raw.loc[raw.warmup.eq(0) & raw.case.str.startswith("syn-")].copy()
    blocks = ["case", "block"]
    if not data.groupby("case")[["input_sha256", "distance"]].nunique().eq(1).all().all():
        raise ValueError("Une configuration doit conserver la même paire et le même score.")
    if not data.groupby(blocks).host.nunique().eq(1).all():
        raise ValueError("Les alignements d’un bloc doivent utiliser la même paire sur le même hôte.")
    times = data.pivot(index=blocks, columns="algorithm", values="seconds")
    if times[["nvpd-cuda", "edlib", "wfa2"]].isna().any().any():
        raise ValueError("Il manque une méthode dans un bloc apparié.")
    ratios = times[["edlib", "wfa2"]].rdiv(times["nvpd-cuda"], axis=0).reset_index()
    parts = ratios.case.str.split("-", expand=True)
    ratios["length"] = parts[1].astype(int)
    ratios["rate"] = (100 * parts[2].astype(float)).round().astype(int)
    return ratios


def _number(value):
    digits = 3 if value < 1 else 2 if value < 10 else 1
    return f"{value:.{digits}f}".replace(".", ",")


def make_figure(records):
    """Présenter les rapports et leurs intervalles sur toute la grille mesurée."""
    use_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    frame = pd.DataFrame(records).set_index(["library", "length", "rate"])
    lengths = sorted(set(frame.index.get_level_values("length")))
    rates = sorted(set(frame.index.get_level_values("rate")))
    colors = {"cpu": COLORS[0], "gpu": COLORS[2], "overlap": COLORS[4]}
    labels = {"cpu": "Bibliothèque CPU plus rapide", "gpu": "NvPD CUDA plus rapide",
              "overlap": "IC recouvrant l’égalité"}
    observed = set()
    groups = (("edlib", "Edlib"), ("wfa2", "WFA2-lib"))
    for j, (library, name) in enumerate(groups):
        for k, rate in enumerate(rates):
            y = j * (len(rates) + .24) + k
            ax.text(-.08, y, f"{rate} %", ha="right", va="center", fontsize=9)
            for x, length in enumerate(lengths):
                row = frame.loc[(library, length, rate)]
                side = "cpu" if row.low95 > 1 else "gpu" if row.high95 < 1 else "overlap"
                observed.add(side)
                ax.add_patch(Rectangle((x + .015, y - .48), .97, .96,
                                       color=colors[side], alpha=.32, lw=0))
                spread = f"{row.sd:.2g}".replace(".", ",")
                ax.text(x + .5, y - .17, f"{_number(row['mean'])} ± {spread}", ha="center",
                        va="center", fontsize=9.6, fontweight="bold")
                ax.text(x + .5, y + .20, f"[{_number(row.low95)} ; {_number(row.high95)}]",
                        ha="center", va="center", fontsize=8.8)
        ax.text(-.60, j * (len(rates) + .24) + (len(rates)-1)/2, name,
                ha="right", va="center", fontweight="bold")
    for i, length in enumerate(lengths):
        unit = f"{length//1000} knt" if length < 1000000 else f"{length//1000000} Mnt"
        ax.text(i+.5, -.76, unit, ha="center", va="center", fontweight="bold")
    ax.text(-.08, -.76, "Taux", ha="right", va="center", fontweight="bold")
    ax.set(xlim=(-1.35, len(lengths)), ylim=(len(rates)*2-.22, -1.))
    ax.set_axis_off()
    fig.legend(handles=[Patch(facecolor=colors[s], edgecolor=INK, alpha=.42, label=labels[s])
                        for s in colors if s in observed], loc="outside upper center", ncol=2)
    return fig
