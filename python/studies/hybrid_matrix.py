"""Temps de la matrice complète et répartition des paires entre Grace et Hopper."""
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, INK, PAPER, column_titles, estimate, label_panels, legend_above, powder_cmap, use_style


COHORTS = {"near": "Un locus", "mixed": "Deux loci", "distant": "Huit loci"}
POLICIES = {
    "cpu-edlib": "Edlib",
    "cpu-wfa2": "WFA2",
    "gpu": "NvPD sur GPU",
    "cpu-sketch": "WFA2 ou Edlib",
    "static": "Partage fixe",
    "hybrid": "Routage strict",
    "hybrid-share": "Routage avec reprise",
}
ROUTED = ("static", "hybrid", "hybrid-share")


def read_trials(path):
    """Conserver les blocs mesurés sur les mêmes entrées pour les sept politiques."""
    data = pd.read_csv(path)
    keys = ["cohort", "block", "policy"]
    if data.duplicated(keys).any() or set(data.cohort) != set(COHORTS) or set(data.policy) != set(POLICIES):
        raise ValueError("Les trois cohortes doivent contenir les sept politiques une fois par bloc.")
    numbers = data[["seconds", "sketch_seconds", "cpu_pairs", "gpu_pairs", "records", "threads"]]
    if not np.isfinite(numbers).all().all() or not data.seconds.gt(0).all():
        raise ValueError("Durée ou décompte manquant dans les blocs.")
    if not (data.sketch_seconds.between(0, data.seconds) & data.threads.eq(4)).all():
        raise ValueError("L’estimation doit être comprise dans le temps total sur quatre cœurs.")
    pairs = data.records * (data.records - 1) // 2
    if not (data.gpu_pairs.between(0, pairs) & data.cpu_pairs.add(data.gpu_pairs).eq(pairs)).all():
        raise ValueError("Le partage doit couvrir le triangle supérieur de la matrice.")
    blocks = sorted(data.block.unique())
    if blocks != list(range(1, len(POLICIES) + 1)):
        raise ValueError("Les sept blocs chronométrés doivent couvrir chaque position dans l’ordre de passage.")
    if (data.groupby(["cohort", "policy"]).position.nunique().ne(len(POLICIES)).any()
            or data.groupby("cohort").matrix_sha256.nunique().ne(1).any()
            or data.matrix_sha256.isna().any()):
        raise ValueError("Ordre incomplet ou matrice différente selon la politique.")
    complete = pd.MultiIndex.from_product([COHORTS, blocks, POLICIES], names=keys)
    if len(complete.difference(data.set_index(keys).index)):
        raise ValueError("Il manque une mesure appariée.")
    return data


def make_figure(trials):
    """Comparer les temps absolus et les décisions des politiques mixtes."""
    use_style()
    fig, (elapsed, assigned) = plt.subplots(1, 2, figsize=(8.05, 3.35),
                                          gridspec_kw={"width_ratios": (1.65, 1)})
    elapsed.set_xscale("log")
    offsets = (-.24, 0, .24)
    for offset, marker, color, (cohort, label) in zip(
            offsets, ("o", "s", "^"), COLORS, COHORTS.items()):
        group = trials[trials.cohort.eq(cohort)].groupby("policy")
        for row, policy in enumerate(POLICIES):
            times = group.get_group(policy).seconds
            estimate(elapsed, row + offset, sample(times), color, marker=marker,
                     horizontal=True, label=label if row == 0 else None)
    means = trials.groupby(["policy", "cohort"]).seconds.mean().unstack()
    for row, policy in enumerate(POLICIES):
        elapsed.plot(means.loc[policy, list(COHORTS)], row + np.asarray(offsets),
                     color=INK, lw=.7, alpha=.65, zorder=1.5)
    elapsed.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 3)))
    elapsed.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}".replace(".", ",")))
    elapsed.set_yticks(np.arange(len(POLICIES)), list(POLICIES.values()))
    elapsed.set_ylim(len(POLICIES) - .45, -.55)
    elapsed.grid(False, axis="y")
    elapsed.set_xlabel("Temps de la matrice [s]")
    elapsed.margins(x=.035)

    groups = trials.groupby(["policy", "cohort"]).gpu_pairs
    counts = groups.mean().unstack().loc[list(ROUTED), list(COHORTS)].to_numpy()
    total = int((trials.cpu_pairs + trials.gpu_pairs).max())
    assigned.pcolormesh(counts, vmin=0, vmax=total, cmap=powder_cmap(),
                        edgecolors=INK, linewidth=.5)
    assigned.set_ylim(len(ROUTED), 0)
    assigned.set_xticks(np.arange(len(COHORTS)) + .5, ["1", "2", "8"])
    assigned.set_yticks(np.arange(len(ROUTED)) + .5, [POLICIES[key] for key in ROUTED])
    assigned.set_xlabel("Nombre de loci")
    assigned.tick_params(length=0)
    assigned.grid(False)
    for row, policy in enumerate(ROUTED):
        for col, cohort in enumerate(COHORTS):
            observed = groups.get_group((policy, cohort))
            avg = counts[row, col]
            if observed.nunique() == 1:
                label = str(int(avg))
            else:
                stats = sample(observed)
                label = (f"{avg:.1f} ± {stats['sd']:.1f}\n"
                         f"[{stats['low95']:.1f} ; {stats['high95']:.1f}]").replace(".", ",")
            assigned.text(col+.5, row+.5, label, ha="center", va="center", fontsize=8,
                          color=PAPER if avg / total > .65 else INK)
    label_panels((elapsed, assigned))
    column_titles((elapsed, assigned), ("Durée totale", "Paires calculées sur GPU"))
    legend_above(fig, (elapsed,), 3)
    return fig


COMPOSITION = {
    "cpu-edlib": "Edlib",
    "cpu-wfa2": "WFA2",
    "hybrid": "Routage par mots",
    "hybrid-counts": "Routage avec effectifs",
}


def read_composition(path):
    """Vérifier les quatre politiques appariées sur les mêmes 64 entrées."""
    data = pd.read_csv(path)
    keys = ["block", "policy"]
    expected = pd.MultiIndex.from_product([range(1, 9), COMPOSITION], names=keys)
    if data.duplicated(keys).any() or not data.set_index(keys).index.sort_values().equals(expected.sort_values()):
        raise ValueError("Chaque politique doit figurer dans les huit blocs appariés.")
    if (not data.records.eq(64).all() or not (data.cpu_pairs + data.gpu_pairs).eq(2016).all()
            or data.matrix_sha256.nunique(dropna=False) != 1 or data.matrix_sha256.isna().any()
            or not np.isfinite(data.seconds).all() or not data.seconds.gt(0).all()
            or not data.sketch_seconds.between(0, data.seconds).all() or not data.threads.eq(4).all()):
        raise ValueError("La comparaison exige la même matrice de 2 016 scores sur quatre cœurs.")
    return data


def composition_figure(data):
    """Montrer le temps de calcul de la matrice des 64 génomes."""
    use_style()
    fig, axis = plt.subplots(figsize=(6.6, 2.05))
    groups = data.groupby("policy").seconds
    for row, key in enumerate(COMPOSITION):
        stats = sample(groups.get_group(key))
        axis.barh(row, stats["mean"], height=.47, color=COLORS[0 if key.startswith("cpu") else 2],
                  edgecolor=INK, linewidth=.6)
        estimate(axis, row, stats, INK, marker="|", horizontal=True)
        label = f"{stats['mean']:.2f}".replace(".", ",")
        axis.annotate(label, (max(stats["high95"], stats["mean"] + stats["sd"]), row),
                      xytext=(5, 0), textcoords="offset points", ha="left", va="center", fontsize=8.5)
    axis.set_yticks(range(len(COMPOSITION)), list(COMPOSITION.values()))
    axis.set_ylim(len(COMPOSITION) - .45, -.55)
    axis.grid(False, axis="y")
    axis.set_xlim(0, axis.get_xlim()[1] * 1.11)
    axis.set_xlabel("Temps de la matrice [s]")
    axis.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}".replace(".", ",")))
    return fig
