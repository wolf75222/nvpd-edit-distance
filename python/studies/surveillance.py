"""Structure de la matrice de surveillance et durée du calcul sur ROMEO."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ct_plot


def ambiguity_bounds(matrix, metadata, content):
    """Borne par différence des effectifs de N pour chaque paire non ordonnée."""
    ids = matrix.index
    if not ids.equals(matrix.columns) or not ids.is_unique:
        raise ValueError("matrix labels do not match")
    groups = metadata.set_index("id").loc[ids, "group"].to_numpy()
    counts = content.set_index("id").loc[ids, "n"].to_numpy()
    i, j = np.triu_indices(len(ids), 1)
    score = matrix.to_numpy()[i, j]
    bound = abs(counts[i] - counts[j])
    if (score <= 0).any() or (bound > score).any():
        raise ValueError("positive distances satisfying the count bound are required")
    subset = np.where((groups[i] == "reference") | (groups[j] == "reference"),
                      "reference", np.where(groups[i] == groups[j], "within", "between"))
    return pd.DataFrame({"group": subset, "bound": bound, "score": score, "ratio": bound / score})


def ambiguity_figure(pairs):
    """Fonctions de survie empiriques, avec les égalités incluses au seuil."""
    ct_plot.use_style()
    fig, ax = plt.subplots(figsize=(6.3, 2.7))
    styles = [("between", "Entre lignages", 0, "-"), ("within", "Même lignage", 1, "--"),
              ("reference", "Avec Wuhan-Hu-1", 2, "-.")]
    for key, name, color, style in styles:
        ratio = np.sort(pairs.loc[pairs.group == key, "ratio"].to_numpy())
        thresholds = np.unique(np.r_[0., ratio, 1.])
        above = 100 * (1 - np.searchsorted(ratio, thresholds) / len(ratio))
        ax.step(thresholds, above, where="pre", color=ct_plot.COLORS[color],
                linestyle=style, lw=2, label=f"{name} ({len(ratio):,})".replace(",", " "))
        ax.plot(.5, 100 * np.mean(ratio >= .5), "o", color=ct_plot.COLORS[color],
                markeredgecolor=ct_plot.INK, markeredgewidth=.5, ms=5)
    ax.axvline(.5, **ct_plot.GUIDE)
    ax.set(xlim=(0, 1), ylim=(0, 101), xticks=np.linspace(0, 1, 6),
           xlabel=r"Seuil $q$ de la borne relative [1]",
           ylabel=r"Paires avec $|\Delta\nu|/d\geq q$ [%]")
    ax.set_xticks(np.linspace(0, 1, 6), ["0", "0,2", "0,4", "0,6", "0,8", "1"])
    ct_plot.legend_above(fig, [ax], ncol=3)
    return fig


def timing_samples(results):
    """Temps des deux allocations de base et comparaison appariée des répartitions."""
    samples = []
    for job, host in [("839288", "epyc"), ("839289", "grace"), ("839304", "epyc")]:
        frame = pd.read_csv(results / f"sars-romeo/{job}/trials.csv").query("block > 0").copy()
        frame["host"] = host
        frame["schedule"] = np.where(frame.threads == 1, "serial", "groups")
        if job == "839288":
            original_split = (frame.algorithm == "nvpd-split") & (frame.threads == 32)
            frame.loc[original_split, "schedule"] = "unfiltered"
            frame = frame[original_split | frame.algorithm.isin(["nvpd-omp", "nvpd-adaptive"])]
        elif job == "839304":
            # Keep both 32-core library variants from their paired A/B experiment.
            frame = frame[(frame.threads == 1) | (frame.algorithm == "nvpd-split")]
        samples.append(frame)
    control = pd.read_csv(results / "sars-schedule/839293/trials.csv").rename(
        columns={"layout": "schedule", "elapsed_s": "wall_s"})
    data = pd.concat([*samples, control.assign(host="epyc")], ignore_index=True)
    grouped = data.groupby(["host", "algorithm", "schedule", "threads"])
    if (data.matrix_sha256.nunique() != 1 or len(grouped) != 18
            or not grouped.block.count().eq(5).all()
            or not grouped.block.nunique().eq(5).all()):
        raise ValueError("expected five complete measured blocks per configuration")
    return data


def timing_figure(records):
    """Temps complets de la matrice, avec le contrôle de répartition sur EPYC."""
    ct_plot.use_style()
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6.5, 4.7),
                             gridspec_kw={"height_ratios": [6, 3]})
    designs = [("epyc", ["edlib", "biwfa", "wfa2", "nvpd-split", "nvpd-adaptive", "nvpd-omp"]),
               ("grace", ["edlib", "wfa2", "nvpd-cuda-local"])]
    styles = {"serial": ("o", 0, "1 cœur"), "rows": ("^", 3, "32 cœurs, rangées"),
              "groups": ("s", 1, "32 cœurs, groupes"), "gpu": ("D", 2, "GH200"),
              "unfiltered": ("X", 5, "32 cœurs, sans filtre")}
    names = {"edlib": "Edlib", "wfa2": "WFA2", "biwfa": "BiWFA",
             "nvpd-split": "Hybride NvPD/BiWFA", "nvpd-adaptive": "NvPD adaptatif",
             "nvpd-omp": "NvPD dense", "nvpd-cuda-local": "NvPD à halo local"}
    for ax, (host, algorithms) in zip(axes, designs):
        for y, algorithm in enumerate(algorithms):
            times = [r for r in records if r["host"] == host and r["algorithm"] == algorithm]
            means = [r["mean"] for r in times]
            if len(means) > 1:
                ax.plot([min(means), max(means)], [y, y], color=ct_plot.COLORS[4], lw=1.25)
            for result in times:
                key = "gpu" if algorithm == "nvpd-cuda-local" else result["schedule"]
                marker, color, label = styles[key]
                ct_plot.estimate(ax, y, result, ct_plot.COLORS[color], marker=marker,
                                 horizontal=True, label=label)
                if host == "grace":
                    left = key == "gpu"
                    ax.annotate(f"{result['mean']:.2f} s".replace(".", ","), (result["mean"], y),
                                xytext=(-8 if left else 8, 0), textcoords="offset points",
                                ha="right" if left else "left", va="center", fontsize=8)
        ax.set_yticks(range(len(algorithms)), [names[k] for k in algorithms])
        ax.set_ylim(len(algorithms)-.55, -.55)
        ax.set_xscale("log")
        ax.grid(False, axis="y")
    axes[-1].set_xlabel("Temps de la matrice [s]")
    ct_plot.label_panels(axes)
    ct_plot.column_titles(axes, ["EPYC 9654", "Grace et Hopper"])
    ct_plot.legend_above(fig, axes, ncol=3)
    return fig
