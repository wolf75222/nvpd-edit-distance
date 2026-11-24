"""Distribution des paires de fenêtres humaines entre rangs MPI."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from ct_plot import COLORS, category_rows, column_titles, estimate, label_panels, legend_above, use_style


def make_figure(root: Path):
    """Tracer l'effet de l'ordre du FASTA à travail identique sur les deux nœuds EPYC."""
    data = root / "results/mpi/schedules/839194"

    def load(file):
        with (data / file).open() as stream:
            return list(csv.DictReader(stream))

    stats = load("summary.csv")
    scale = load("scaling.csv")
    policies = ("static", "cyclic", "hilbert", "hilbert-weighted", "steal")
    labels = ("Blocs", "Cyclique", "Hilbert", "Hilbert pondéré", "Vol")
    metrics = (("gcups", "Débit du lot [GCUPS]", 1), ("max_mean_rank_cells", "Travail maximal / moyen [1]", 1),
               ("busy_fraction", "Temps occupé des fils [%]", 100), ("efficiency", "Efficacité forte [%]", 100))
    layouts = (("locus", "Ordre des loci", "o", -.15), ("length", "Ordre des longueurs", "s", .15))

    use_style()
    fig, axes = plt.subplots(2,2,figsize=(8.15,5.1),layout="constrained")
    for axis,(metric,xlabel,factor) in zip(axes.flat,metrics):
        table = scale if metric == "efficiency" else [r for r in stats if r["nodes"] == "2"]
        lookup = {(r["layout"],r["schedule"]):r for r in table if r["metric"] == metric}
        for i,policy in enumerate(policies):
            for j,(layout,legend,marker,offset) in enumerate(layouts):
                row = lookup[layout,policy]
                center = {k:float(row[k])*factor for k in ("mean","sd","low95","high95")}
                estimate(axis,i+offset,center,COLORS[j],marker=marker,label=legend,horizontal=True)
        category_rows(axis,labels,xlabel)
        low, high = axis.get_xlim()
        if low < 0 or (metric != "max_mean_rank_cells" and low < .12*high):
            axis.set_xlim(left=0)
    column_titles(axes[0],("Deux nœuds", "Deux nœuds"))
    column_titles(axes[1],("Deux nœuds", "Un à deux nœuds"))
    label_panels(axes.flat)
    legend_above(fig,axes.flat,2)
    return fig
