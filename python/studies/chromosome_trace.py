"""Coût de la reconstruction du CIGAR et gain dû au découpage des chromosomes."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from ct_plot import COLORS, INK, PAPER, category_rows, column_titles, estimate, label_panels, legend_above, use_style


def _summary(path):
    with path.open(newline="") as stream:
        return {(row["case"], row["method"], int(row["threads"]), row["metric"]):
                {key: float(row[key]) for key in ("mean", "sd", "high95")}
                for row in csv.DictReader(stream) if row["mode"] == "trace"}


def make_figure(root: Path):
    base = root / "results/trace/chromosomes"
    whole = _summary(base / "839200/summary.csv")
    ablation = _summary(base / "839193/summary.csv")
    use_style()
    fig, (cost, gain) = plt.subplots(1, 2, figsize=(9.1, 2.85), width_ratios=(1, 1.28))
    methods = [("split", 8, "Fragments · 8 cœurs", (6, 1)), ("wfa2", 1, "WFA2", (6, -3)),
               ("biwfa", 1, "BiWFA", (6, -1)), ("edlib", 1, "edlib", (-6, -10))]
    for i, (method, threads, label, offset) in enumerate(methods):
        t = whole["chr1", method, threads, "seconds"]
        m = whole["chr1", method, threads, "peak_rss_kib"]
        x, y = t["mean"], m["mean"] / 1024
        for key, width, cap in (("sd", 3.0, 0), ("high95", .8, 2)):
            dx = t[key] if key == "sd" else t[key] - x
            dy = (m[key] if key == "sd" else m[key] - m["mean"]) / 1024
            cost.errorbar(x, y, xerr=dx, yerr=dy, fmt="none", ecolor=COLORS[i], elinewidth=width, capsize=cap)
        cost.plot(x, y, marker="os^D"[i], ls="none", color=COLORS[i], markeredgecolor=INK)
        cost.annotate(label, (x, y), xytext=offset, textcoords="offset points", fontsize=8,
                      va="center", ha="right" if offset[0] < 0 else "left",
                      bbox=dict(facecolor=PAPER, edgecolor="none", pad=0.8))
    cost.set(xscale="log", xlabel="durée de reconstruction [s]", ylabel="pic résident [Mio]",
             xlim=(0.006, 300), ylim=(300, 2550), yticks=[500, 1000, 1500, 2000, 2500])
    cases = ("chr21-2m", "chr20", "chr18", "chr1")
    metrics = [("speedup_trim-only", "paire entière · 1 cœur", "s"),
               ("speedup_split", "fragments · 1 cœur", "o")]
    for i, case in enumerate(cases):
        for j, (metric, label, marker) in enumerate(metrics):
            stat = ablation[case, "split", 8, metric]
            row = i + (j - 0.5) * .22
            gain.hlines(row, 1, stat["mean"], color=COLORS[j+1], linewidth=2.2)
            estimate(gain, row, stat, COLORS[j+1], marker=marker, label=label, horizontal=True)
    gain.axvline(1, color=INK, linewidth=.7)
    gain.set(xscale="log", xlim=(0.8, 54))
    gain.set_xticks([1, 2, 4, 8, 16, 32], ["1", "2", "4", "8", "16", "32"])
    category_rows(gain, ["Fenêtre chr21", "Chr20", "Chr18", "Chr1"], "temps de référence / temps à huit cœurs [1]")
    column_titles((cost, gain), ("CIGAR du chromosome 1", "Effet du découpage et des cœurs"))
    label_panels((cost, gain))
    legend_above(fig, [gain], ncol=2)
    return fig
