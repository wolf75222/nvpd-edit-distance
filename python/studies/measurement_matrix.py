"""Matrice de mesures répétées. La teinte reprend une couleur de ct_plot."""

import numpy as np
from matplotlib.patches import Rectangle

from ct_plot import COLORS, INK


def _number(value, precision):
    if 0 < abs(value) < .001:
        return "<0,001"
    if precision.endswith("g"):
        # A fixed decimal representation keeps large timings readable on a printed page.
        digits = int(precision[1:-1])
        places = digits - 1 - int(np.floor(np.log10(abs(value)))) if value else 0
        value = round(value, places)
        return f"{value:,.{max(0, places)}f}".replace(",", " ").replace(".", ",")
    if precision.endswith("f") and 0 < value < .5 * 10 ** (-int(precision[1:-1])):
        return "<" + format(10 ** (-int(precision[1:-1])), precision).replace(".", ",")
    return format(value, precision).replace(".", ",")


def measurement_matrix(axis, records, rows, columns, *, precision=".2f", spread=".2g",
                       color=0, log=False, center=None, reverse=False):
    """Inscrire la moyenne, l'écart-type et la demi-largeur de l'IC dans chaque case."""
    transform = np.log10 if log else float
    means = [transform(entry["mean"]) for entry in records.values()]
    low, high = min(means), max(means)
    axis.grid(False)
    axis.set(xlim=(-.5, len(columns) - .5), ylim=(len(rows) - .5, -.5),
             xticks=range(len(columns)), xticklabels=columns,
             yticks=range(len(rows)), yticklabels=rows)
    axis.tick_params(length=0, pad=4)
    for row in range(len(rows)):
        for col in range(len(columns)):
            entry = records.get((row, col))
            if entry is None:
                axis.text(col, row, "n.m.", ha="center", va="center", fontsize=8)
                continue
            mean = entry["mean"]
            value = transform(mean)
            if center is None:
                fraction = (value - low) / (high - low) if high > low else 0
                if reverse:
                    fraction = 1 - fraction
                shade = COLORS[color]
            else:
                extent = max(abs(high - center), abs(low - center))
                fraction = abs(value - center) / extent if extent else 0
                shade = COLORS[2 if value >= center else 1]
            axis.add_patch(Rectangle((col-.495, row-.495), .99, .99,
                                     facecolor=shade, alpha=.10 + .63*fraction, linewidth=0))
            axis.text(col, row-.16, _number(mean, precision), color=INK, ha="center",
                      va="center", fontweight="bold", fontsize=8.3)
            detail = " · ".join(_number(v, spread) for v in
                                (entry["sd"], entry["high95"] - mean))
            axis.text(col, row+.19, detail, ha="center", va="center", fontsize=6.9)
