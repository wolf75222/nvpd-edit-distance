"""Tableaux compacts pour les mesures dont les unités diffèrent selon les colonnes."""

import numpy as np

from ct_plot import INK


def booktabs(axis, headers, rows, *, widths):
    """Aligner les libellés à gauche et centrer les valeurs dans les colonnes."""
    if not (len(headers) == len(widths) and all(len(row) == len(headers) for row in rows)):
        raise ValueError("Chaque ligne doit avoir une entrée par colonne.")
    bounds = np.r_[0, np.cumsum(widths)] / sum(widths)
    count = len(rows)
    axis.set(xlim=(0, 1), ylim=(-.05, count + 1.18))
    axis.set_axis_off()
    for i, row in enumerate([headers, *rows]):
        y = count + .63 if i == 0 else count - i + .49
        for j, value in enumerate(row):
            x = bounds[j] + .014 if j == 0 else (bounds[j] + bounds[j + 1]) / 2
            axis.text(x, y, value, fontsize=10, va="center", ha="left" if j == 0 else "center",
                      color=INK, fontweight="bold" if i == 0 else "normal", linespacing=1.12)
    axis.hlines([count + 1.15, count + .09, -.005], 0, 1, color=INK, linewidths=[1, .5, 1])
