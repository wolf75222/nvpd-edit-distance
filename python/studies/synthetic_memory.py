"""Empreinte des 60 paires de la grille et RSS échantillonné à paire fixée."""
import matplotlib.pyplot as plt
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, INK, column_titles, estimate, label_panels, legend_above, use_style


LENGTHS = (5000, 20000, 100000, 1000000)
LENGTH_LABELS = ("5 knt", "20 knt", "100 knt", "1 Mnt")
RATES = (2, 5, 10)


def read_memory(root):
    directory = root / "results"
    campaign = pd.read_csv(directory / "romeo/summary.csv")
    gpu = campaign.loc[campaign.algorithm.eq("nvpd-cuda")].copy()
    gpu["percent"] = (100 * gpu.rate).round().astype(int)
    grid = pd.MultiIndex.from_product([LENGTHS, RATES, range(1, 6)])
    if (not gpu.set_index(["length", "percent", "rep"]).index.sort_values().equals(grid)
            or gpu.input_sha256.nunique() != 60):
        raise ValueError("La grille attend cinq FASTA distincts par longueur et taux.")
    expected = 24 * (gpu[["m", "n"]].max(axis=1) + 1)
    if not gpu.max_device_bytes.eq(expected).all():
        raise ValueError("Le volume GPU diffère des deux lignes et des quatre tables d'occurrence.")
    gpu["device_mib"] = gpu.max_device_bytes / 2**20
    gpu["host_mib"] = gpu.max_peak_rss_kib / 1024

    raw = pd.read_csv(directory / "repeated/gh200/memory_trials.csv")
    libraries = raw[raw.case.str.startswith("syn-") & raw.algorithm.isin(["wfa2", "edlib"])]
    libraries = libraries.merge(gpu[gpu.rep.eq(1)][["input_sha256", "length", "percent"]],
                                on="input_sha256", how="left", validate="many_to_one")
    observations = libraries.groupby(["length", "percent", "algorithm"]).block
    expected = pd.Series([9 if length <= 20000 else 5
                          for length, _, _ in observations.size().index], index=observations.size().index)
    if (len(expected) != 24 or libraries.length.isna().any()
            or not libraries.rss_metric.eq("sampled_smaps_rollup_max").all()
            or not observations.size().eq(expected).all() or not observations.nunique().eq(expected).all()):
        raise ValueError("La grille de RSS requiert un maximum échantillonné par processus.")
    libraries["rss_mib"] = libraries.rss_kib / 1024
    return gpu, libraries


def make_figure(gpu, libraries):
    use_style()
    fig = plt.figure(figsize=(7.5, 3.85))
    grid = fig.add_gridspec(2, 3, width_ratios=(1.45, 1, 1))
    footprint = fig.add_subplot(grid[:, 0])
    axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in (1, 2)]
    shapes = gpu.groupby("length")
    for field, label, color, marker in (
            ("host_mib", "RSS hôte NvPD", COLORS[4], "s"),
            ("device_mib", "Tampons CUDA", COLORS[2], "D")):
        stats = shapes[field].agg(["mean", "min", "max"]).loc[list(LENGTHS)]
        values = stats["mean"]
        footprint.plot(range(4), values, color=color, marker=marker, markeredgecolor=INK,
                       label=label)
        footprint.vlines(range(4), stats["min"], stats["max"], color=INK)
        for pos, value in enumerate(values):
            if field == "host_mib" and pos in (1, 2):
                continue
            digits = 0 if field == "host_mib" else (1 if value > 10 else 2)
            text = f"{value:.{digits}f}".replace(".", ",")
            footprint.annotate(text, (pos, value), xytext=(0, 7), textcoords="offset points",
                               ha="center", fontsize=8)
    footprint.set_yscale("log")
    footprint.set(ylim=(.08, 430), xlim=(-.27, 3.32), xticks=range(4),
                  xticklabels=LENGTH_LABELS, ylabel="mémoire [Mio]", xlabel="longueur de référence")
    footprint.set_yticks([.1, 1, 10, 100], ["0,1", "1", "10", "100"])
    footprint.grid(False, axis="x")

    for axis, length in zip(axes, LENGTHS):
        selected = libraries[libraries.length.eq(length)]
        for row, rate in enumerate(RATES):
            values = selected[selected.percent.eq(rate)].groupby("algorithm").rss_mib
            points = {algorithm: sample(counts) for algorithm, counts in values}
            axis.hlines(row, points["edlib"]["mean"], points["wfa2"]["mean"],
                        color=COLORS[4], lw=1)
            for algorithm, label, color, marker in (("edlib", "edlib", COLORS[0], "o"),
                                                     ("wfa2", "WFA2-lib", COLORS[1], "s")):
                estimate(axis, row, points[algorithm], color, marker=marker, label=label, horizontal=True)
        axis.set(yticks=range(3), yticklabels=["2 %", "5 %", "10 %"], ylim=(2.48, -.47))
        axis.locator_params(axis="x", nbins=4)
        axis.grid(False, axis="y")
        axis.margins(x=.13)
    for axis in axes[2:]:
        axis.set_xlabel("RSS sondé [Mio]")
    label_panels([footprint, *axes])
    column_titles([footprint, *axes], ["NvPD sur GH200", *LENGTH_LABELS])
    legend_above(fig, [footprint, *axes], ncol=4)
    return fig
