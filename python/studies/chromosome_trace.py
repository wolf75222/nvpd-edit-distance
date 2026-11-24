"""Temps, mémoire et accélération du calcul par fragments sur quatre paires humaines."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from ct_plot import COLORS, category_rows, column_titles, estimate, label_panels, legend_above, use_style


def make_figure(root: Path):
    dataset = root / "results/trace/chromosomes/839187/summary.csv"
    with dataset.open() as stream:
        stats = {(r["case"],r["mode"],r["method"],int(r["threads"]),r["metric"]):r
                 for r in csv.DictReader(stream)}
    with (root/"results/trace/chromosomes/839200/summary.csv").open() as stream:
        for row in csv.DictReader(stream):
            if row["metric"] in ("seconds","peak_rss_kib"):
                stats[row["case"],row["mode"],row["method"],int(row["threads"]),row["metric"]] = row
    cases = ("chr21-2m","chr20","chr18","chr1")
    names = ("Fenêtre chr21","Chr20","Chr18","Chr1")
    methods = (("split",8,"Fragments (8 cœurs)"),("wfa2",1,"WFA2"),
               ("biwfa",1,"BiWFA"),("edlib",1,"edlib"))

    def draw(axis, position, row, scale, index, label):
        if row:
            values = {k:float(row[k])*scale for k in ("mean","sd","high95")}
            estimate(axis,position,values,COLORS[index % len(COLORS)],marker=("o","s","D","^","v","P","X")[index],
                     label=label,horizontal=True)

    use_style()
    fig, axes = plt.subplots(2,2,figsize=(8.35,4.9),layout="constrained")
    for ax, mode, metric, scale, unit in (
            (axes[0,0],"score","seconds",1000,"Temps [ms]"),
            (axes[0,1],"trace","seconds",1000,"Temps [ms]"),
            (axes[1,0],"trace","peak_rss_kib",1/1024,"Pic résident [Mio]")):
        for i,case in enumerate(cases):
            for j,(backend,threads,label) in enumerate(methods):
                draw(ax,i+(j-1.5)*.14,stats.get((case,mode,backend,threads,metric)),scale,j,label)
        ax.set_xscale("log")
        category_rows(ax,names,unit)
    ax = axes[1,1]
    for i,case in enumerate(cases):
        for j,threads in enumerate((2,4,8)):
            draw(ax,i+(j-1)*.16,stats[case,"trace","split",threads,"speedup_split"],
                 1,j+4,f"{threads} cœurs")
    category_rows(ax,names,"Accélération par rapport à un cœur [1]")
    column_titles(axes[0],("Score global","CIGAR complet"))
    column_titles(axes[1],("Mémoire du CIGAR","Parallélisme des fragments"))
    label_panels(axes.flat)
    legend_above(fig,axes.flat,4)
    return fig
