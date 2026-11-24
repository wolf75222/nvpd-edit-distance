"""Rapports de temps par lot et transferts vers le CPU avec WFA-GPU."""
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

from ct_plot import COLORS, INK, column_titles, estimate, join, label_panels, legend_above, use_style


def _observations(path):
    times = defaultdict(dict)
    fallbacks = defaultdict(set)
    with path.open() as stream:
        for row in csv.DictReader(stream):
            block = int(row['block'])
            if block == 0:
                continue
            if int(row['wrong_scores']):
                raise ValueError(f"{row['case']}: scores discordants")
            key = (row['case'], row['method'], int(row['cap']))
            if block in times[key]:
                raise ValueError(f'{key}: bloc en double')
            times[key][block] = float(row['seconds'])
            fallbacks[key].add(int(row['cpu_fallback']))
    return times, fallbacks


def _ratio(numerator, denominator):
    if numerator.keys() != denominator.keys() or set(numerator) != set(range(1, 6)):
        raise ValueError('La comparaison requiert cinq blocs complets et appariés.')
    values = [t / denominator[i] for i, t in numerator.items()]
    mean, sd = statistics.mean(values), statistics.stdev(values)
    return dict(mean=mean, sd=sd, high95=mean + 2.7764451052 * sd / 5**.5)


def make_figure(root):
    """Comparer le temps de processus WFA-GPU à WFA2 sur un cœur Grace."""
    use_style()
    data = Path(root) / 'results/wfa-gpu'
    synthetic, phi_fallbacks = _observations(data / '839291/trials.csv')
    _, human_fallbacks = _observations(data / '839294/trials.csv')
    fig, (ax, lower) = plt.subplots(2, 1, figsize=(7.0, 4.1), height_ratios=[1, 1.15])
    ax.set_xscale('log')
    ax.set_xlim(.55, 140)
    ax.set_xticks([.6, 1, 2, 5, 10, 25, 50, 100], ['0,6', '1', '2', '5', '10', '25', '50', '100'])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1, color=INK, lw=1)
    labels, positions = [], []
    for y, (length, rate, cap) in enumerate(((5000,'0.02',108), (5000,'0.1',508),
                                               (20000,'0.02',408), (20000,'0.1',2008))):
        labels.append(f'{length//1000} knt  ·  {round(100*float(rate))} %')
        positions.append(y)
        centers = []
        for size, marker, color in ((5,'o',COLORS[0]), (128,'s',COLORS[1])):
            case = f'n{length}_r{rate}_p{size}'
            result = _ratio(synthetic[case,'wfa-gpu',cap], synthetic[case,'wfa2',cap])
            centers.append(result['mean'])
            estimate(ax, y, result, color, marker=marker, horizontal=True,
                     label=f'{size} paires' if y == 0 else None)
        join(ax, [y, y], centers, COLORS[4], horizontal=True)
    ax.set_yticks(positions, labels)
    ax.set_ylim(3.5, -.5)
    ax.set_xlabel('Rapport des temps de processus WFA-GPU / WFA2-lib [1]')

    samples = [('φX174', 'phix174_all', 256, 2016, phi_fallbacks)]
    samples.extend((f'Chr. 21, K = {k:,}'.replace('1,024','1\u202f024'), 'chr21', k, 768, human_fallbacks)
                   for k in (64,256,1024))
    allocations = []
    for name, case, cap, count, counts in samples:
        key = (case, 'wfa-gpu', cap)
        if len(counts[key]) != 1:
            raise ValueError(f'{case}: nombre de reprises variable entre répétitions')
        fallbacks, = counts[key]
        allocations.append((name, fallbacks, count))
    with (data / '839298/trials.csv').open() as stream:
        matrix = [row for row in csv.DictReader(stream) if row['method']=='hybrid' and int(row['block'])>0]
    if len(matrix) != 5 or any(int(row['wrong_scores']) for row in matrix):
        raise ValueError('La matrice SARS-CoV-2 requiert cinq blocs concordants.')
    routes = {(int(row['cpu_scored']), int(row['pairs'])) for row in matrix}
    if len(routes) != 1:
        raise ValueError('La répartition des paires SARS-CoV-2 change entre répétitions.')
    allocations.insert(1, ('SARS-CoV-2', *routes.pop()))
    for i, (_, cpu, count) in enumerate(allocations):
        percentage = 100 * cpu / count
        lower.barh(i, percentage, height=.57, color=COLORS[2], edgecolor=INK, linewidth=.6,
                   zorder=2, label='CPU' if i == 0 else None)
        lower.barh(i, 100-percentage, left=percentage, height=.57, color=COLORS[3], edgecolor=INK,
                   linewidth=.6, alpha=.52, zorder=2, label='GPU' if i == 0 else None)
        inside = percentage > 75
        lower.annotate(f'{cpu:,}/{count:,}'.replace(',', '\u202f'), (percentage, i),
                       xytext=(-4 if inside else 4, 0), ha='right' if inside else 'left',
                       textcoords='offset points', va='center', fontsize=8)
    lower.grid(False, axis='y')
    lower.set_yticks(range(len(allocations)), [item[0] for item in allocations])
    lower.set_ylim(len(allocations)-.5, -.5)
    lower.set_xlim(0,100)
    lower.set_xticks([0,25,50,75,100])
    lower.set_xlabel('Répartition du calcul [% des paires]')
    label_panels([ax,lower])
    column_titles([ax,lower], ['Coût relatif sur les lots synthétiques', 'Répartition des paires biologiques'])
    legend_above(fig, [ax,lower], 4)
    return fig
