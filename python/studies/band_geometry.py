"""Cellules visitées et chemin d'édition d'un intervalle du chromosome 20."""
from collections import Counter
import csv
import re

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

from ct_plot import COLORS, INK, PAPER, column_titles, label_panels, legend_above, use_style


def _trace(cigar, left, right):
    i = j = 0
    path, edits = [(i, j)], []
    for count, op in re.findall(r"(\d+)([=XID])", cigar):
        for _ in range(int(count)):
            if op in "=X" and (left[i] == right[j]) != (op == "="):
                raise ValueError("Le CIGAR contredit les deux séquences.")
            i += op != "I"
            j += op != "D"
            if op != "=":
                edits.append(((path[-1][0] + i) / 2, (path[-1][1] + j) / 2))
            path.append((i, j))
    if (i, j) != (len(left), len(right)):
        raise ValueError("Le CIGAR doit consommer les deux fragments.")
    return np.array(path), np.array(edits)


def study(root):
    folder = root / "results/adaptive/geometry"
    _, left, _, right = (folder / "chr20.fa").read_text().splitlines()
    with (folder / "check.csv").open() as stream:
        rows = {row["algorithm"]: row for row in csv.DictReader(stream)}
    work = rows["nvpd-adaptive"]
    prefix, suffix = int(work["prefix"]), int(work["suffix"])
    a, b = (x[prefix:len(x)-suffix] for x in (left, right))
    m, n = len(a), len(b)
    total = int(work["cells"])
    first = max(1, sum((Counter(b) - Counter(a)).values()))
    last = int(work["threshold"])
    if int(work["passes"]) != 2 or m != n:
        raise ValueError("Cette reconstruction requiert deux passages et deux longueurs égales.")
    i, j = np.mgrid[1:m+1, 1:n+1]
    narrow, wide = (abs(j-i) <= k//2 for k in (first, last))
    # La seconde passe arrive à la dernière ligne. Le reste du compteur
    # donne le nombre cumulé de mises à jour jusqu'à l'arrêt de la première.
    cutoff, = np.flatnonzero(narrow.sum(axis=1).cumsum() == total - wide.sum()) + 1
    visits = wide.astype(int) + (narrow & (i <= cutoff)).astype(int)
    path, edits = _trace(rows["edlib"]["cigar"], left, right)
    score = int(work["distance"])
    if len(edits) != score or score != int(rows["edlib"]["distance"]):
        raise ValueError("Le coût du chemin et les deux scores doivent coïncider.")
    path = path - prefix
    path = path[((path >= 0) & (path <= [m, n])).all(axis=1)]
    return dict(visits=visits, path=path, edits=edits-prefix, first=first, last=last,
                cutoff=int(cutoff))


def make_figure(data):
    use_style()
    visits, path = data["visits"], data["path"]
    m, n = visits.shape
    k, stop = data["first"], data["cutoff"]
    margin = data["last"] // 2
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 2.65), gridspec_kw={"width_ratios": [1, 2.3]})
    locus, band = axes
    cmap = ListedColormap([PAPER, COLORS[0], COLORS[1]])
    locus.imshow(visits.T, origin="lower", extent=(.5,m+.5,.5,n+.5), vmin=0, vmax=2,
                 cmap=cmap, interpolation="nearest")
    locus.plot(*path.T, color=INK, lw=.7)
    locus.set(xlabel="Préfixe de A [nt]", ylabel="Préfixe de B [nt]",
              xlim=(-.5,m+.5), ylim=(-.5,n+.5), xticks=[0,9,18,27], yticks=[0,9,18,27])
    locus.text(.07, .94, f"{m*n} cases", transform=locus.transAxes, va="top")

    offsets = np.arange(-margin, margin+1)
    lines = np.arange(1, m+1)
    strip = np.zeros((len(offsets), m), dtype=int)
    for row, offset in enumerate(offsets):
        columns = lines+offset
        valid = (1 <= columns) & (columns <= n)
        strip[row, valid] = visits[lines[valid]-1, columns[valid]-1]
    band.pcolormesh(np.arange(.5, m+1.5), np.arange(-margin-.5, margin+1.5), strip,
                    cmap=cmap, vmin=0, vmax=2, edgecolors=PAPER, linewidth=.5, antialiased=True)
    band.plot(path[:, 0], path[:, 1]-path[:, 0], color=INK, lw=1.35, label="Chemin optimal")
    edits = data["edits"]
    band.plot(edits[:, 0], edits[:, 1]-edits[:, 0], ls="none", marker="o", ms=4,
              markerfacecolor=PAPER, markeredgecolor=INK, label="Édition")
    for value, label in [(1, "Une évaluation"), (2, "Deux évaluations")]:
        band.plot([], [], lw=0, marker="s", markersize=7, markeredgewidth=0,
                  markerfacecolor=cmap(value/2), label=label)
    band.set(xlabel="Longueur du préfixe de A, $i$ [nt]", ylabel="Décalage $j-i$ [nt]",
             xticks=[0,5,10,15,20,25,27], yticks=offsets, xlim=(-.5,m+.5), ylim=(-3.1,3.6))
    band.annotate(f"Arrêt à $i={stop}$ pour $k={k}$", (stop, k//2+.45),
                  xytext=(stop-1.5, 2.75), ha="center", va="bottom",
                  arrowprops=dict(arrowstyle="->", color=INK, lw=.75))
    for ax in axes:
        ax.grid(False)
    label_panels(axes)
    column_titles(axes, ("Grille élaguée", "Bande dépliée par diagonale"))
    legend_above(fig, [band], 4)
    return fig
