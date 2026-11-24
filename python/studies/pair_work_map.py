"""Carte du travail et des coupes de Hilbert sur le panel du chromosome 21."""
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.path import Path
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

from ct_plot import COLORS, INK, PAPER, column_titles, label_panels, legend_above, powder_cmap, use_style


def _hilbert_order(side, i, j):
    """Même orientation du carré et mêmes rotations que dans scheduled_matrix.hpp."""
    x, y = i.copy(), j.copy()
    keys = np.zeros_like(x)
    scale = side // 2
    while scale:
        rx, ry = (x & scale) != 0, (y & scale) != 0
        keys += scale**2 * ((3 * rx) ^ ry)
        rotate = ~ry
        reflect = rx & rotate
        x[reflect], y[reflect] = side - 1 - x[reflect], side - 1 - y[reflect]
        x[rotate], y[rotate] = y[rotate], x[rotate]
        scale //= 2
    return np.argsort(keys)


def work_map(root):
    data = root / "results/mpi/schedules"
    fragments = pd.read_csv(data / "windows.tsv", sep="\t").sort_values(["length", "id"])
    lengths = fragments.length.to_numpy(dtype=np.int64)
    n, ranks = len(lengths), 48
    if n != 128 or fragments.id.nunique() != n or not np.all(fragments.end - fragments.start == lengths):
        raise ValueError("Les identifiants ou les coordonnées des 128 fragments sont incohérents.")
    i, j = np.triu_indices(n, k=1)
    order = _hilbert_order(n, i, j)
    i, j = i[order], j[order]
    work = lengths[i] * lengths[j]
    prefix = np.r_[0, np.cumsum(work)]
    plain = len(work) * np.arange(ranks + 1) // ranks
    targets = int(prefix[-1]) * np.arange(ranks + 1) // ranks
    upper = np.searchsorted(prefix, targets)
    lower = np.maximum(0, upper - 1)
    weighted = np.where(targets - prefix[lower] < prefix[upper] - targets, lower, upper)

    raw = pd.read_csv(data / "839194/rank_work.csv")
    trials = raw[raw.nodes.eq(2) & raw.layout.eq("length")]
    owners = np.full((n, n), -1, dtype=int)
    top = (i >= 3*n//4) & (j >= 3*n//4)
    metrics = {
        "pairs": int(top.sum()), "pairs_total": len(work), "cells": int(work[top].sum()),
        "cells_total": int(prefix[-1]), "min_length": int(lengths[3*n//4]), "max_length": int(lengths[-1])
    }
    for policy, cuts in (("hilbert", plain), ("hilbert-weighted", weighted)):
        count, cells = np.diff(cuts), np.diff(prefix[cuts])
        group = trials[trials.schedule.eq(policy)]
        if set(group.block) != set(range(1, 6)):
            raise ValueError(f"Blocs incomplets pour {policy}.")
        for _, block in group.groupby("block"):
            observed = block.sort_values("rank")
            if (not np.array_equal(observed["rank"], np.arange(ranks))
                or not np.array_equal(observed.pairs, count)
                or not np.array_equal(observed.cells, cells)):
                raise ValueError(f"Les coupes diffèrent des comptes par rang pour {policy}.")
        owner = np.repeat(np.arange(ranks), count)
        if policy == "hilbert":
            owners[j, i] = owner
        else:
            owners[i, j] = owner
        metrics[policy] = int(np.unique(owner[top]).size)
    costs = np.multiply.outer(lengths, lengths) / 1e6
    np.fill_diagonal(costs, np.nan)
    return costs, owners, metrics


def _boundaries(owners):
    segments = []
    for direction, grid in enumerate((owners, owners.T)):
        mask = (grid[:, :-1] != grid[:, 1:]) & (grid[:, :-1] >= 0) & (grid[:, 1:] >= 0)
        y, x = np.where(mask)
        edges = np.stack((np.column_stack((x + 1.5, y + .5)), np.column_stack((x + 1.5, y + 1.5))), axis=1)
        segments.append(edges if direction == 0 else edges[:, :, ::-1])
    return np.concatenate(segments)


def make_figure(costs, owners):
    use_style()
    n = len(costs)
    first = 3*n//4 + 1
    fig, axes = plt.subplots(1, 2, figsize=(6.35, 3.80), layout="constrained")
    edges = _boundaries(owners)
    image = axes[0].imshow(costs, origin="lower", interpolation="nearest", cmap=powder_cmap(),
                           norm=LogNorm(np.nanmin(costs), np.nanmax(costs)), alpha=.52,
                           extent=(.5, n+.5, .5, n+.5))
    zoom = owners[first-1:, first-1:]
    y, x = np.indices(zoom.shape)
    upper = y > x
    classes = np.ma.masked_where(x == y, upper.astype(int))
    axes[1].imshow(classes, origin="lower", interpolation="nearest",
                   cmap=ListedColormap([COLORS[2], COLORS[1]]), vmin=0, vmax=1, alpha=.26,
                   extent=(first-.5, n+.5, first-.5, n+.5))
    for triangle in (upper, y < x):
        for rank in np.unique(zoom[triangle]):
            region = triangle & (zoom == rank)
            distance = distance_transform_edt(np.pad(region, 1))[1:-1, 1:-1]
            candidates = np.argwhere(distance == distance.max())
            center = np.argwhere(region).mean(axis=0)
            row, column = candidates[np.argmin(np.sum((candidates - center)**2, axis=1))]
            axes[1].text(column+first, row+first, str(rank), ha="center", va="center", fontsize=8.6)
    for axis, start, ticks in zip(axes, (1, first), ([1, 32, 64, 96, 128], [97, 104, 112, 120, 128])):
        axis.add_collection(LineCollection(edges, colors=INK, linewidths=.4, alpha=.84))
        axis.plot([.5, n+.5], [.5, n+.5], color=PAPER, lw=1.1)
        axis.set(xlim=(start-.5, n+.5), ylim=(start-.5, n+.5), xticks=ticks, yticks=ticks,
                 xlabel="Indice du fragment [1]", aspect="equal")
        axis.grid(False)
    axes[0].set_ylabel("Indice du fragment [1]")
    axes[0].add_patch(Rectangle((first-.5, first-.5), n-first+1, n-first+1,
                                facecolor="none", edgecolor=COLORS[1], linewidth=1.6, zorder=4))
    for vertices, color, text in (([(-1, -1), (-1, 1), (1, 1)], COLORS[1], "Haut · même effectif cible"),
                                  ([(-1, -1), (1, -1), (1, 1)], COLORS[2], "Bas · même charge cible")):
        axes[0].plot([], [], marker=Path([*vertices, vertices[0]]), ls="none",
                     markerfacecolor=color, markeredgecolor=INK, markersize=8, label=text)
    legend_above(fig, axes, 2)
    label_panels(axes)
    column_titles(axes, ("Travail sur le panel entier", "Rangs des paires longues"))
    scale = fig.colorbar(image, ax=axes, location="bottom", shrink=.77, aspect=42, pad=.04,
                         ticks=[.3, 1, 4, 16, 64, 150])
    scale.ax.set_xticklabels(["0,3", "1", "4", "16", "64", "150"])
    scale.ax.xaxis.set_minor_locator(NullLocator())
    scale.set_label("Travail par paire dans le panneau a) [millions de cellules]")
    return fig
