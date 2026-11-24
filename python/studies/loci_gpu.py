"""Temps du score global pour des fenêtres homologues et des loci distincts de chr21."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, INK, estimate, legend_above, use_style

JOB = '839323'
CASES = ('matched-262144', 'different-131072', 'different-262144')
METHODS = ('wfa2', 'edlib', 'nvpd-cuda-local')


def read_trials(root: Path):
    table = pd.read_csv(root / 'results' / 'energy' / 'gh200' / JOB / 'loci.csv')
    data = table[table.block > 0].copy()
    cells = data.groupby(['case', 'algorithm']).block
    if (len(data) != 45 or cells.ngroups != 9 or not cells.nunique().eq(5).all()
            or set(data.block) != set(range(1, 6)) or set(data.case) != set(CASES)
            or set(data.algorithm) != set(METHODS) or not data.seconds.gt(0).all()
            or not table.groupby('case').distance.nunique().eq(1).all()):
        raise ValueError('incomplete chromosome candidate comparison')
    return data


def make_figure(root: Path):
    data = read_trials(root)
    use_style()
    figure, axis = plt.subplots(figsize=(7.4, 2.8), layout='constrained')
    region_labels = ['Homologue · 262 144 nt', 'Loci distincts · 131 072 nt',
                     'Loci distincts · 262 144 nt']
    speeds = []
    for row, case in enumerate(CASES):
        subset = data[data.case == case]
        stats = {key: sample(subset.loc[subset.algorithm == key, 'seconds']) for key in METHODS}
        means = [result['mean'] for result in stats.values()]
        axis.hlines(row, min(means), max(means), color=INK, linewidth=.65, zorder=1.4)
        for name, label, color, marker in zip(METHODS, ('WFA2', 'Edlib', 'NvPD compacté'),
                                               COLORS, ('s', 'D', 'o')):
            result = stats[name]
            estimate(axis, row, result, color, marker=marker, label=label if not row else None,
                     horizontal=True)
            value = result['mean']
            time = (f'{value*1000:.2f} ms' if value < .01 else f'{value:.2f} s').replace('.', ',')
            axis.annotate(time, (value, row), xytext=(0, 6), textcoords='offset points',
                          ha='center', va='bottom', fontsize=8)
            speeds.extend([result['low95'], result['high95']])
    axis.set(xscale='log', xlim=(min(speeds) / 1.9, max(speeds) * 1.8),
             ylim=(2.40, -.55), yticks=np.arange(3), yticklabels=region_labels,
             xlabel="Temps de calcul du score [s]")
    axis.minorticks_off()
    axis.grid(visible=False, axis='y')
    legend_above(figure, (axis,), 3)
    return figure
