"""Traduction d’adresses et réemploi des chromosomes CHM13 sur EPYC."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, GUIDE, INK, category_rows, column_titles, estimate, label_panels, legend_above, use_style

JOB = "839316"


def load_pages(root: Path):
    path = root / "results/profile/huge-pages" / JOB
    trials = pd.read_csv(path / "trials.csv")
    trials["walks"] = trials.dtlb_walks / trials.calls
    trials["score_ms"] = 1000 * trials.core_s / trials.calls
    reuse = pd.read_csv(path / "reuse.csv").pivot(index="block", columns=["queries", "pages"], values="total_s")
    delta = (reuse.xs("base", level="pages", axis=1) - reuse.xs("thp", level="pages", axis=1)) * 1000
    return trials, delta


def make_figure(root: Path):
    trials, delta = load_pages(root)
    use_style()
    fig, (walk, cost) = plt.subplots(2, 1, figsize=(6.9, 4.45), height_ratios=[1.05, 1])
    labels = []
    for solver, name in (("split", "Fragments"), ("wfa2", "WFA2")):
        for case in ("chr1", "chr18", "chr20"):
            case_data = trials[(trials.case == case) & (trials.solver == solver)]
            stats = {p: sample(g.walks) for p, g in case_data.groupby("pages")}
            row = len(labels)
            walk.hlines(row, stats["thp"]["mean"], stats["base"]["mean"], color=COLORS[4], linewidth=1.1)
            for p, key, color, marker in (("base", "4 Kio", COLORS[0], "s"), ("thp", "2 Mio", COLORS[1], "o")):
                estimate(walk, row, stats[p], color, marker=marker, label=key, horizontal=True)
            labels.append(f"{name} · {case}")
    category_rows(walk, labels, "Parcours de tables demandés par score [1]")
    walk.set_xscale("log")
    walk.set(xlim=(1, 220000))
    walk.set_xticks([1, 10, 100, 1000, 10000, 100000], ["1", "10", "100", "1 000", "10 000", "100 000"])

    split = trials[(trials.case == "chr1") & (trials.solver == "split")]
    paired = split.pivot(index="block", columns="pages", values=["setup_s", "score_ms"])
    offset = (paired.setup_s.base - paired.setup_s.thp) * 1000
    step = paired.score_ms.base - paired.score_ms.thp
    queries = np.linspace(1, delta.columns.max(), 80)
    model = [sample(offset + k*step) for k in queries]
    cost.fill_between(queries, [r["low95"] for r in model], [r["high95"] for r in model],
                      color=COLORS[2], alpha=.4, lw=0)
    cost.plot(queries, [r["mean"] for r in model], color=INK, lw=1.2, label="Modèle linéaire")
    cost.axhline(0, **GUIDE)
    for k, times in delta.items():
        estimate(cost, k, sample(times), COLORS[5], marker="D", label="Série chronométrée")
    cost.set(xlabel="Scores successifs sur les mêmes entrées [1]",
             ylabel="Temps à 4 Kio moins\nle temps à 2 Mio [ms]", xlim=(-8, delta.columns.max() + 10),
             xticks=[1, 96, 192, 288, 384])
    column_titles((walk, cost), ("Traduction des adresses", "Coût cumulé du solveur par fragments sur le chromosome 1"))
    label_panels((walk, cost))
    legend_above(fig, [walk, cost], ncol=4)
    return fig
