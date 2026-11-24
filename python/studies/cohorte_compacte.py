"""Distance des haplotypes du chromosome 21 et temps des matrices sur ROMEO."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import t as student

import ct_plot
from bench_stats import sample as summarize

METHODS = {
    "edlib": "Edlib",
    "wfa2": "WFA2",
    "guarded": "Deux bits, bords testés",
    "one-sided": "Deux bits, un front",
    "byte": "Octets, deux fronts",
    "packed": "Deux bits, deux fronts",
}


def read_cohort(root):
    data = Path(root) / "results/human/packed-cohort"
    def read(pattern):
        frames = [pd.read_csv(p) for p in sorted(data.glob(pattern))]
        if not frames:
            raise ValueError(f"Aucun relevé pour {pattern}")
        return pd.concat(frames, ignore_index=True).query("trial > 0").copy()
    trials, weak = read("trials-*.csv"), read("weak-[0-9]*.csv")
    distribution = pd.read_csv(data / "weak-distance-profile.csv").set_index("threads")
    weak["mean_squared_distance"] = weak.threads.map(distribution.mean_squared_distance)
    trials["mean_squared_distance"] = distribution.at[64, "mean_squared_distance"]
    hist = pd.read_csv(data / "distances.csv")
    if hist.pairs.sum() != 1024*1023//2:
        raise ValueError("Distribution de distances incomplète.")
    for frame in (trials, weak):
        for _, group in frame.groupby(["threads", "method", "collapse"]):
            if sorted(group.trial) != list(range(1, 6)):
                raise ValueError("La série doit contenir cinq blocs.")
    if set(trials.threads) != {1, 2, 4, 8, 16, 32, 64} or set(weak.threads) != {1, 2, 4, 8, 16, 32}:
        raise ValueError("Plan de mesure incomplet.")
    return trials, weak, hist


def efficiency_estimates(trials, weak):
    fixed = trials.query("method == 'packed' and collapse == 0")
    proportional = pd.concat([weak, fixed.query("threads == 64")])
    for label, frame in [("Forte", fixed), ("Faible", proportional), ("Pondérée", proportional)]:
        reference = frame.query("threads == 1")
        base = summarize(reference.total_s)
        for threads, group in frame.groupby("threads"):
            if threads == 1:
                yield label, threads, dict(mean=100., sd=0., high95=100.)
                continue
            measured = summarize(group.total_s)
            work_ratio = 1. if label == "Forte" else group.kernel_calls.iloc[0] / reference.kernel_calls.iloc[0]
            if label == "Pondérée":
                work_ratio *= group.mean_squared_distance.iloc[0] / reference.mean_squared_distance.iloc[0]
            ratio = 100 * work_ratio * base["mean"] / (threads * measured["mean"])
            terms = [(s["sd"]/s["mean"])**2 for s in (base, measured)]
            sd = ratio * sum(terms)**.5
            # Five blocks per configuration; rank counts run in separate processes.
            df = 4 * sum(terms)**2 / sum(v*v for v in terms)
            yield label, threads, dict(mean=ratio, sd=sd, high95=ratio + student.ppf(.975, df)*sd / 5**.5)


def make_figure(trials, weak, hist):
    ct_plot.use_style()
    fig = plt.figure(figsize=(7.05, 5.55))
    layout = fig.add_gridspec(2, 2, height_ratios=[1.23, 1], width_ratios=[1, 1])
    times = fig.add_subplot(layout[0, :])
    density = fig.add_subplot(layout[1, 0])
    scaling = fig.add_subplot(layout[1, 1])
    for y, (method, label) in enumerate(METHODS.items()):
        centers = []
        for collapse, color, marker, name in [(False, ct_plot.COLORS[0], "o", "Matrice directe"),
                                               (True, ct_plot.COLORS[1], "s", "Séquences regroupées")]:
            cases = trials.query("threads == 64 and method == @method and collapse == @collapse")
            res = summarize(cases.total_s * 1000)
            centers.append(res["mean"])
            ct_plot.estimate(times, y, res, color, marker=marker, horizontal=True,
                             label=name if not y else None)
        ct_plot.join(times, [y,y], centers, ct_plot.COLORS[4], horizontal=True)
    times.set_xscale("log")
    times.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000],
                     ["10", "20", "50", "100", "200", "500", "1 000", "2 000", "5 000"])
    times.set_yticks(range(len(METHODS)), list(METHODS.values()))
    times.set_ylim(len(METHODS)-0.5, -0.5)
    times.set_xlabel("Temps de la matrice [ms]")
    density.bar(hist.distance, 100*hist.pairs/hist.pairs.sum(), width=1,
                color=ct_plot.COLORS[0], zorder=3)
    density.set_xlim(-1, 101)
    density.set_ylim(bottom=0)
    density.set_xlabel("Distance globale [éditions]")
    density.set_ylabel("Paires [%]")
    series = list(efficiency_estimates(trials, weak))
    for label, color, marker, key in [("Forte", ct_plot.COLORS[2], "^", "Forte"),
                                      ("Faible", ct_plot.COLORS[5], "D", "Faible par paire"),
                                      ("Pondérée", ct_plot.COLORS[3], "v", r"Faible par $d^2$")]:
        groups = [(p, values) for name, p, values in series if name == label]
        scaling.plot([p for p,_ in groups], [est["mean"] for _,est in groups], color=color, lw=1)
        for p, est in groups:
            ct_plot.estimate(scaling, p, est, color, marker=marker, label=key if p == 1 else None)
    scaling.axhline(100, **ct_plot.GUIDE)
    scaling.set_xscale("log", base=2)
    ticks = [1, 2, 4, 8, 16, 32, 64]
    scaling.set_xticks(ticks, [str(p) for p in ticks])
    scaling.set_xlabel("Cœurs physiques")
    scaling.set_ylabel("Efficacité [%]")
    ct_plot.column_titles([times, density, scaling],
                          ["Coût sur un nœud", "Diversité locale", "Changement d’échelle"])
    ct_plot.label_panels([times, density, scaling])
    ct_plot.legend_above(fig, [times, scaling], ncol=5)
    return fig
