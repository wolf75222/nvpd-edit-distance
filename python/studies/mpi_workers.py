"""Durées de la matrice humaine avec les deux moteurs de travailleurs MPI."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from ct_plot import COLORS, GUIDE, category_rows, column_titles, estimate, label_panels, legend_above, use_style


def make_figure(root: Path):
    with (root/"results/mpi/workers/839197/summary.csv").open() as stream:
        data = {(r["schedule"],int(r["nodes"]),r["backend"],r["metric"]):r for r in csv.DictReader(stream)}
    use_style()
    fig,axes = plt.subplots(1,2,figsize=(8.25,2.85),layout="constrained")
    policies = ("static","cyclic","hilbert","hilbert-weighted","steal")
    labels = ("Blocs","Cyclique","Hilbert","Hilbert pondéré","Vol")
    for i,policy in enumerate(policies):
        series = ((0,1,"kokkos","speedup_over_openmp","192 cœurs",0),
                  (0,2,"kokkos","speedup_over_openmp","384 cœurs",1),
                  (1,2,"openmp","seconds","OpenMP",2), (1,2,"kokkos","seconds","Kokkos",3))
        for panel,nodes,backend,metric,name,index in series:
            row = data[policy,nodes,backend,metric]
            stat = {k:float(row[k]) for k in ("mean","sd","high95")}
            estimate(axes[panel],i+(-.13 if index%2==0 else .13),stat,COLORS[index],
                     marker=("o","s","D","^")[index],label=name,horizontal=True)
    category_rows(axes[0],labels,"Rapport de temps OpenMP / Kokkos [1]")
    category_rows(axes[1],labels,"Temps [s]")
    axes[0].axvline(1,**GUIDE)
    column_titles(axes,("Effet du moteur de travailleurs","Matrice sur deux nœuds"))
    label_panels(axes)
    legend_above(fig,axes,4)
    return fig
