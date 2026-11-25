"""Comparaison appariée des récurrences de NvPD sur les cœurs EPYC."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, GUIDE, INK, PAPER, column_titles, label_panels, legend_above, use_style


def block_times(trials):
    """Aligner les mesures par bloc et vérifier les scores, les temps et les effectifs."""
    keys = ["cohort", "length", "threads", "block"]
    if trials.duplicated(keys + ["variant"]).any():
        raise ValueError("Plusieurs mesures pour le même bloc et la même variante.")
    if (not np.isfinite(trials[["seconds", "warmup_seconds"]]).all().all()
            or not trials[["seconds", "warmup_seconds"]].gt(0).all().all()
            or trials.groupby(["cohort", "length"]).score.nunique().ne(1).any()):
        raise ValueError("Temps ou score invalide dans la comparaison des récurrences.")
    values = trials.pivot(index=keys, columns="variant", values="seconds")
    if values.empty or values.isna().any().any() or len(values.columns) < 2:
        raise ValueError("Chaque bloc doit contenir les mêmes variantes.")
    nblocks = values.groupby(level=["cohort", "length", "threads"]).size()
    if nblocks.nunique() != 1 or nblocks.min() < 3:
        raise ValueError("Nombre de blocs incomplet pour une des configurations.")
    return values


def paired_gains(trials, reference):
    """Former les rapports de temps dans chaque bloc, à longueur et affinité fixées."""
    values = block_times(trials)
    return values.drop(columns=reference).rdiv(values[reference], axis="index").reset_index()


def validation_times(trials, length, reference, selected):
    """Conserver les temps mesurés sur une même longueur et calculer leur rapport par bloc."""
    times = block_times(trials[trials.length.eq(length)])[[reference, selected]].copy()
    times.columns = ["initial", "integrated"]
    times["gain"] = times.initial / times.integrated
    return times.reset_index()


VARIANTS = {
    "wide": "SIMD 512 bits",
    "centered": "Scores centrés",
    "centered-wide": "Centrage et SIMD 512 bits",
    "blocked": "Permutations par blocs",
}


FOLLOW = {
    "index32": "Indices de 32 bits",
    "aligned": "Chargements alignés",
    "unroll2": "Déroulage ×2",
    "unroll4": "Déroulage ×4",
}


def block_logs(gains, variants):
    """Un logarithme de gain par variante et par bloc, en pondérant les deux paires à parts égales."""
    keys = ["length", "threads", "block"]
    if gains.groupby(keys).cohort.nunique().ne(2).any():
        raise ValueError("Chaque bloc de sélection doit apparier les deux fenêtres.")
    logs = gains[list(variants)].apply(np.log)
    return pd.concat([gains[keys], logs], axis=1).groupby(keys, as_index=False).mean()


def selection_logs(gains):
    return block_logs(gains[gains.length.eq(32768)], VARIANTS)


def geometric_interval(axis, y, values, color, marker, label):
    """Reporter l’écart-type et l’IC calculés dans l’espace des logarithmes."""
    stats = sample(values)
    mean = np.exp(stats["mean"])
    axis.hlines(y, 1, mean, color=color, lw=1.3, alpha=.7, zorder=1)
    for width, cap, alpha, half in (
        (3.2, 0, .52, stats["sd"]), (.85, 2.4, 1, stats["high95"] - stats["mean"])
    ):
        left, right = np.exp(stats["mean"] + np.array([-half, half]))
        axis.errorbar(mean, y, xerr=[[mean - left], [right - mean]], fmt="none", ecolor=color,
                      elinewidth=width, capsize=cap, alpha=alpha)
    axis.plot(mean, y, ls="none", color=color, marker=marker, markeredgecolor=INK,
              markeredgewidth=.45, label=label, zorder=3)


def validation_table(axis, times):
    """Publier dans chaque case la moyenne, l’écart-type et la demi-largeur de l’IC."""
    cases = (("near", 1), ("distant", 1), ("near", 4), ("distant", 4))
    axis.set_xlim(-.50, 3.56)
    axis.set_ylim(3.55, -1.30)
    axis.set_axis_off()
    for x, name in zip((-.46, 1.21, 2.18, 3.13), ("Paire et cœurs", "Initial [s]", "Intégré [s]", "Gain")):
        axis.text(x, -.86, name, ha="left" if x < 0 else "center", va="center",
                  fontweight="bold", fontsize=8.3)
    for y, width in ((-1.22, .85), (-.55, .55), (1.50, .38), (3.52, .85)):
        axis.hlines(y, -.48, 3.55, color=INK, lw=width)
    for row, (cohort, threads) in enumerate(cases):
        group = times[times.cohort.eq(cohort) & times.threads.eq(threads)]
        locus = "Même locus" if cohort == "near" else "Loci distincts"
        axis.text(-.46, row, f"{locus}\n{threads} {'cœur' if threads == 1 else 'cœurs'}",
                  ha="left", va="center", fontsize=8)
        for x, column in zip((1.21, 2.18, 3.13), ("initial", "integrated", "gain")):
            stat = sample(group[column])
            avg = f"{stat['mean']:.3f}".replace(".", ",")
            spread = f"{stat['sd']:.3f} · {stat['high95'] - stat['mean']:.3f}".replace(".", ",")
            axis.text(x, row-.15, avg, ha="center", va="center", fontweight="bold", fontsize=8.2)
            axis.text(x, row+.18, spread, ha="center", va="center", fontsize=7.3)


def make_figure(logs, times):
    """Distinguer le critère de sélection des temps du binaire sur les séquences prolongées."""
    use_style()
    fig, (training, validation) = plt.subplots(1, 2, figsize=(8.45, 2.8),
                                              gridspec_kw={"width_ratios": (1, 1.15)})
    groups = logs.groupby("threads")
    for threads, dy, color, marker in ((1, -.14, COLORS[0], "o"), (4, .14, COLORS[1], "s")):
        label = "Un cœur" if threads == 1 else "Quatre cœurs"
        for row, variant in enumerate(VARIANTS):
            geometric_interval(training, row+dy, groups.get_group(threads)[variant],
                               color, marker, label if row == 0 else None)
    training.axvline(1, **GUIDE)
    training.set_yticks(np.arange(len(VARIANTS)), list(VARIANTS.values()))
    training.set_ylim(len(VARIANTS) - .46, -.5)
    training.set_xlabel("Temps initial / temps de la variante")
    training.grid(False, axis="y")
    training.margins(x=.06)
    validation_table(validation, times)
    label_panels((training, validation))
    length_label = f"{times.length.unique().item():,}".replace(",", " ")
    column_titles((training, validation),
                  ("Sélection à 32 768 bases", f"Temps à {length_label} bases"))
    legend_above(fig, (training,), 2)
    return fig


def follow_figure(logs):
    """Gains des variantes par rapport au noyau par blocs sur quatre configurations."""
    use_style()
    fig, axis = plt.subplots(figsize=(6.60, 2.88))
    groups = list(logs.groupby(["length", "threads"], sort=True))
    axis.set_xlim(-.5, len(FOLLOW)-.5)
    axis.set_ylim(len(groups)-.5, -.5)
    axis.set_xticks(range(len(FOLLOW)), FOLLOW.values(), fontweight="bold")
    axis.xaxis.tick_top()
    axis.set_yticks(range(len(groups)), [f"{length:,} bases\n{threads} cœur{'s' if threads > 1 else ''}"
                                       .replace(",", " ") for (length, threads), _ in groups])
    axis.grid(False)
    axis.tick_params(axis="both", which="major", length=0, pad=5)
    # Le fond exprime la position de l'intervalle entier par rapport à l'unité.
    colors = (COLORS[1], PAPER, COLORS[2])
    legend = ("IC inférieur à 1", "IC contient 1", "IC supérieur à 1")
    for row, (_, group) in enumerate(groups):
        for col, variant in enumerate(FOLLOW):
            s = sample(group[variant])
            mean, sd, lo, hi = np.exp([s[k] for k in ("mean", "sd", "low95", "high95")])
            category = 0 if hi < 1 else 2 if lo > 1 else 1
            axis.add_patch(Rectangle((col-.5, row-.5), 1, 1, facecolor=colors[category],
                                     alpha=.65, edgecolor=INK, linewidth=.45, label=legend[category]))
            axis.text(col, row-.17, f"{mean:.3f}".replace(".", ","), ha="center", va="center",
                      fontweight="bold", fontsize=9)
            interval = f"[{lo:.3f} ; {hi:.3f}]".replace(".", ",")
            variation = f"$s_g$ = {sd:.3f}".replace(".", ",")
            axis.text(col, row+.16, f"{interval}\n{variation}", ha="center", va="center",
                      linespacing=.99, fontsize=7.5)
    legend_above(fig, (axis,), len(colors))
    return fig
