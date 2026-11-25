"""Temps de deux routages sur quinze paires de confirmation à 100 knt."""
import argparse
import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[4] / "python"))
import ct_plot


GROUPS = ("2pct", "5pct", "10pct")
METHODS = ("hybrid-original", "hybrid-refined", "wfa2")
LABELS = ("Routage initial", "Routage révisé", "WFA2")
COLORS = (ct_plot.COLORS[0], ct_plot.COLORS[2], ct_plot.COLORS[4])


def _seconds(path):
    with path.open() as handle:
        return [row for row in csv.DictReader(handle) if row["metric"] == "seconds"]


def _format(value):
    return f"{value:.2f}".replace(".", ",")


def figure(analysis_dir):
    """Tracer les temps par paire et les IC des moyennes entre cinq paires."""
    directory = Path(analysis_dir)
    pairs = _seconds(directory / "pairs.csv")
    summary = {(r["group"], r["method"]): r for r in _seconds(directory / "summary.csv")}
    ct_plot.use_style()
    graph, (axis, table) = plt.subplots(1, 2, figsize=(6.7, 3.36),
                                        width_ratios=(2.85, 1.15), sharey=True)
    positions = (0, 1.25, 2.5)
    for center, group in zip(positions, GROUPS):
        for index, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS)):
            samples = sorted((int(r["rep"]), 1000 * float(r["mean"])) for r in pairs
                             if r["group"] == group and r["method"] == method)
            if [rep for rep, _ in samples] != list(range(12, 17)):
                raise ValueError(f"Confirmation incomplète : {group}, {method}")
            y = center + (index - 1)*.30
            data = summary[group, method]
            mean, low, high = [1000 * float(data[k]) for k in ("mean", "low95", "high95")]
            axis.barh(y, mean, height=.22, color=color, edgecolor=ct_plot.INK,
                      linewidth=.45, label=label if group == GROUPS[0] else None)
            axis.plot([value for _, value in samples], [y-.085+.022*i for i in range(len(samples))],
                      linestyle="none", marker=".", color=ct_plot.INK, markersize=2.7)
            axis.errorbar(mean, y+.045, xerr=[[mean-low], [high-mean]], color=ct_plot.INK,
                          capsize=2, elinewidth=.85, fmt="none")
            table.text(.24, y, _format(mean), ha="right", va="center",
                       transform=table.get_yaxis_transform())
            table.text(.64, y, f"[{_format(low)} ; {_format(high)}]", ha="center", va="center",
                       transform=table.get_yaxis_transform())
    axis.set(yticks=positions, yticklabels=("2 %", "5 %", "10 %"),
             ylim=(positions[-1]+.53, positions[0]-.53), xlim=(0, 104),
             xlabel="Temps par paire [ms]", ylabel="Éditions introduites")
    axis.yaxis.grid(False)
    table.axis("off")
    for x, text in ((.22, "Moyenne"), (.64, "IC à 95 %")):
        table.annotate(text, (x, 1), xycoords="axes fraction", xytext=(0, 2),
                       textcoords="offset points", ha="center", va="bottom", fontweight="bold")
    ct_plot.legend_above(graph, (axis,), ncol=3)
    return graph


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    result = figure(arguments.analysis)
    result.savefig(arguments.output)
    plt.close(result)
