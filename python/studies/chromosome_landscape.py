"""Localisation des révisions de CHM13 dans les coordonnées de la version 1.0."""
from pathlib import Path

import numpy as np
import pandas as pd


def extract_intervals(root, vcf):
    """Extraire les intervalles d'allèles du VCF du consortium, dont le rejeu contrôle l'empreinte."""
    from runpy import run_path

    panel = root / "results/human/chromosome-panel"
    names = pd.read_csv(panel / "manifest.tsv", sep="\t").case
    replay = run_path(str(root / "bench/human/chromosome-panel/release_patch.py"))
    variants = replay["read_variants"](vcf, set(names))
    rows = []
    for chrom in names:
        for start, ref, alt in variants[chrom]:
            rows.append((chrom, start, start + len(ref), replay["allele_edit_cost"](ref, alt)))
    return pd.DataFrame(rows, columns=["case", "start0", "end0", "allele_cost"])


def read_panel(root):
    panel = root / "results/human/chromosome-panel"
    lengths = pd.read_csv(panel / "manifest.tsv", sep="\t").set_index("case").m
    records = pd.read_csv(panel / "release-intervals.csv").join(lengths, on="case", validate="many_to_one")
    if not ((0 <= records.start0) & (records.start0 < records.end0) & (records.end0 <= records.m)).all():
        raise ValueError("Intervalle hors du chromosome de référence.")
    records["right_width"] = records.m - records.start0
    records["radius"] = np.minimum(records.end0, records.right_width)
    archived = pd.read_csv(panel / "release-patch.csv").set_index("case")
    for observed, count, cost in [(records, "records", "path_cost"),
                                  (records[records.radius <= 1000], "terminal_records", "terminal_cost")]:
        grouped = observed.groupby("case").allele_cost
        if not (grouped.size().reindex(lengths.index, fill_value=0).eq(archived[count]).all()
                and grouped.sum().reindex(lengths.index, fill_value=0).eq(archived[cost]).all()):
            raise ValueError("Les intervalles ne reproduisent pas le bilan du rejeu VCF.")
    return records, lengths


def region_costs(records, chromosomes, width):
    side = np.select([records.end0 <= width, records.right_width <= width],
                     ["left", "right"], default="other")
    costs = records.groupby([side, "case"]).allele_cost.sum().unstack()
    return costs.reindex(index=["left", "right", "other"], columns=chromosomes).fillna(0).astype(int)


def make_figure(records, costs, width):
    import matplotlib.pyplot as plt
    from ct_plot import COLORS, INK, PAPER, column_titles, label_panels, legend_above, powder_cmap, use_style

    use_style()
    fig, (profile, grid) = plt.subplots(2, 1, figsize=(7.4, 4.25), height_ratios=(2.6, 1))
    by_radius = records.groupby("radius").allele_cost.agg(["size", "sum"])
    for key, label, color, style in [("size", "Entrées du VCF", COLORS[0], "--"),
                                      ("sum", "Coût des allèles", COLORS[1], "-")]:
        series = by_radius[key].cumsum() / by_radius[key].sum() * 100
        crop = series.loc[:1000]
        profile.step(np.r_[0, crop.index, 1000], np.r_[0, crop.values, crop.iloc[-1]],
                     where="post", color=color, ls=style, lw=1.9, label=label)
        share = crop.loc[:width].iloc[-1]
        profile.plot(width, share, "o", color=color, markeredgecolor=INK, zorder=3)
        label_value = f"{share:.1f} %".replace(".", ",")
        profile.annotate(label_value, (width, share), xytext=(-6, 1) if key == "sum" else (8, -11),
                         ha="right" if key == "sum" else "left",
                         textcoords="offset points", bbox=dict(facecolor=PAPER, edgecolor="none", pad=1.0))
    profile.axvline(width, ls=":", color=INK, lw=.8)
    largest = records.loc[records.allele_cost.idxmax()]
    radius = largest.radius
    before = by_radius.loc[by_radius.index < radius, "sum"].sum()
    middle = 100 * (before + by_radius.loc[radius, "sum"]/2) / by_radius["sum"].sum()
    profile.annotate(f"{largest['case']} · {largest.allele_cost} éditions",
                     xy=(radius, middle), xytext=(.20, .835), textcoords="axes fraction",
                     ha="left", va="center", fontsize=9,
                     arrowprops=dict(arrowstyle="->", color=INK, lw=.75, shrinkB=2),
                     bbox=dict(facecolor=PAPER, edgecolor="none", pad=1.5))
    profile.set(xlim=(0, 1000), ylim=(0, 105), xticks=[0, 250, 500, 750, 1000],
                yticks=[0, 25, 50, 75, 100],
                xlabel="Largeur de la marge à chaque extrémité [nt]", ylabel="Part cumulée [%]")
    legend_above(fig, [profile], 2)
    columns = len(costs.columns)
    max_log_cost = np.log1p(costs.to_numpy().max())
    grid.grid(False)
    grid.pcolormesh(np.arange(columns+1), np.arange(4), np.log1p(costs.to_numpy()),
                    cmap=powder_cmap(), vmin=0, vmax=2*max_log_cost, edgecolors=PAPER, linewidth=.6)
    for i, j in np.ndindex(costs.shape):
        value = costs.iat[i, j]
        if value:
            grid.text(j+.5, i+.5, str(value), ha="center", va="center", fontsize=8.5)
    grid.set(ylim=(3, 0), yticks=[.5, 1.5, 2.5], yticklabels=["Début", "Fin", "Autres"],
             xticks=np.arange(columns)+.5, xticklabels=[c.removeprefix("chr") for c in costs.columns],
             xlabel="Chromosome de CHM13 v1.0")
    grid.tick_params(length=0, pad=2.5)
    column_titles([profile, grid], ["Couverture des révisions depuis les extrémités",
                                     f"Répartition du coût pour des marges de {width} nt [éditions]"])
    label_panels([profile, grid])
    return fig


if __name__ == "__main__":
    import argparse

    cli = argparse.ArgumentParser(description="Extraire les intervalles du VCF de révision de CHM13.")
    cli.add_argument("vcf", type=Path)
    cli.add_argument("output", type=Path)
    args = cli.parse_args()
    extract_intervals(Path(__file__).resolve().parents[2], args.vcf).to_csv(args.output, index=False)
