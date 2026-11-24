"""Travail de la récurrence NvPD, temps et modèle mémoire empirique sur EPYC."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ct_plot import COLORS, GUIDE, INK, column_titles, estimate, label_panels, use_style


def make_figure(root: Path):
    """Comparer les sept entrées et les étalons mesurés du job EPYC."""
    data=root/"results/adaptive/839188"
    def load(file):
        with (data/file).open() as source:
            return list(csv.DictReader(source))
    def uncertainty(row, factor=1):
        return {k:float(row[k])*factor for k in ("mean","sd","low95","high95")}
    runs={(r["case"],r["algorithm"],r["metric"]):r for r in load("summary.csv")}
    calib={(r["kernel"],int(r["footprint_kib"])):r for r in load("calibration-summary.csv")}
    pmu={(r["case"],r["algorithm"],r["metric"]):r for r in load("pmu-summary.csv")}
    work={r["case"]:r for r in load("work.csv")}
    cases=("syn-2","syn-10","independent","chr21-25k","chr21-100k","chr21-250k","chr21-1000k")
    labels=("2 %","10 %","Indép.","25 k","100 k","250 k","1 M")
    use_style()
    fig,axes=plt.subplots(2,2,figsize=(8.2,6.8),layout="constrained")
    time,ratio,bandwidth,roof=axes.flat
    methods=(("nvpd-adaptive","NvPD adaptatif"),("nvpd-omp","NvPD dense"),
             ("banded-trim","Bande classique"),("edlib","edlib"),("wfa2","WFA2"),("biwfa","BiWFA"))
    for j,(alg,title) in enumerate(methods):
        for i,case in enumerate(cases):
            row=runs.get((case,alg,"seconds"))
            if row:
                estimate(time,i+(j-2.5)*.11,uncertainty(row,1000),COLORS[j],
                         marker=("o","s","^","v","D","P")[j],label=title if i==0 else None)
    time.set_yscale("log");time.set_ylabel("Temps du score [ms]")
    time.legend(loc="lower center",bbox_to_anchor=(.5,1.12),ncol=3,fontsize=6)
    for i,case in enumerate(cases):
        ratio.plot(i-.19,float(work[case]["cell_ratio"]),"^",color=COLORS[3],
                   label="Cellules denses / visitées" if i==0 else None)
        for j,(alg,title) in enumerate((("nvpd-omp","Temps dense / adaptatif"),
                                         ("band-control","Temps scalaire / SIMD"))):
            row=runs.get((case,alg,"paired_time_over_adaptive"))
            if row:
                estimate(ratio,i+.02+.16*j,uncertainty(row),COLORS[j],marker=("o","s")[j],
                         label=title if i==0 else None)
    ratio.axhline(1,**GUIDE);ratio.set_yscale("log")
    ratio.set_ylabel("Rapport [1]");ratio.legend(loc="lower center",bbox_to_anchor=(.5,1.12),ncol=2,fontsize=6)
    for ax in (time,ratio):
        ax.set_xticks(range(len(cases)),labels,rotation=22)
        ax.axvline(2.5,color=INK,lw=.55)
        ax.set_xlabel("Trois paires synthétiques, puis quatre fenêtres de chr21")
    sizes=(12,384,12288,196608)
    for j,(kind,text) in enumerate((("read","Lecture"),("triad","Triade"))):
        positions=np.arange(4)+(j-.5)*.10
        for i,(pos,size) in enumerate(zip(positions,sizes)):
            estimate(bandwidth,pos,uncertainty(calib[kind,size]),COLORS[j],marker=("o","s")[j],
                     label=text if i==0 else None)
    bandwidth.set_xticks(range(4),("12 Kio","384 Kio","12 Mio","192 Mio"))
    bandwidth.set_xlabel("Taille cumulée des tableaux actifs");bandwidth.set_ylabel("Débit utile [Go/s]")
    bandwidth.set_ylim(bottom=0);bandwidth.legend(loc="lower center",bbox_to_anchor=(.5,1.12),ncol=2,fontsize=6)
    peak=float(calib["row",12]["mean"]);bw=float(calib["read",196608]["mean"])
    points=[]
    for j,(case,name) in enumerate((("syn-10","Synthétique à 10 %"),("independent","Indépendantes"),
                                     ("chr21-100k","Chr21 à 100 k"))):
        ix=uncertainty(pmu[case,"nvpd-adaptive","updates_per_estimated_dram_byte"])
        iy=uncertainty(pmu[case,"nvpd-adaptive","gcell_s"])
        points.append(ix["mean"])
        roof.errorbar(ix["mean"],iy["mean"],xerr=ix["high95"]-ix["mean"],
                      fmt="none",ecolor=COLORS[j],elinewidth=.85,capsize=2.4)
        estimate(roof,ix["mean"],iy,COLORS[j],marker=("o","s","^")[j],label=name)
    xx=np.geomspace(peak/bw/8,max(points)*2.2,260)
    p_sd=float(calib["row",12]["sd"]);b_sd=float(calib["read",196608]["sd"])
    roof.fill_between(xx,np.minimum(peak-p_sd,(bw-b_sd)*xx),np.minimum(peak+p_sd,(bw+b_sd)*xx),
                      color=INK,alpha=.09,linewidth=0)
    roof.plot(xx,np.minimum(peak,bw*xx),label="Repère étalonné",**GUIDE)
    roof.set_xscale("log");roof.set_yscale("log");roof.set_ylim(peak/12,peak*2)
    roof.set_xlabel("Intensité de lecture estimée [cellules/octet]")
    roof.set_ylabel("Débit NvPD adaptatif [Gcell/s]");roof.legend(loc="lower center",bbox_to_anchor=(.5,1.12),ncol=2,fontsize=6)
    column_titles(axes[0],("Appels complets","Réduction du travail et des temps"))
    column_titles(axes[1],("Débits de lecture et de triade","Modèle de débit en lecture"))
    label_panels(axes.flat)
    return fig
