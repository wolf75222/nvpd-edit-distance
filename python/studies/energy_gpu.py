"""Rendement des tuiles NvPD et coût de lecture des indices sur le chromosome 21."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bench_stats import sample
from ct_plot import COLORS, GUIDE, column_titles, estimate, label_panels, legend_above, use_style
from studies.measurement_matrix import measurement_matrix

DEPTH_JOB = '839310'
PACK_JOB = '839317'
WINDOW_JOB = '839319'


def read_trials(root: Path):
    base = root / 'results' / 'energy' / 'gh200'
    depths = pd.read_csv(base / DEPTH_JOB / 'trials.csv')
    packed = pd.read_csv(base / PACK_JOB / 'packing.csv')
    for frame in (depths, packed):
        frame['time_s'] = frame.seconds / frame.iterations
        frame['energy_j'] = frame.joules / frame.iterations
        frame['gcups_w'] = frame.m * frame.n / (1e9 * frame.energy_j)
    expected = depths.groupby('case').distance.first()
    if (len(packed) != 40 or packed.duplicated(['case', 'block', 'kernel']).any()
            or not packed.distance.eq(packed.case.map(expected)).all()
            or not set(packed.kernel) == {'packed', 'direct'}):
        raise ValueError('inconsistent paired GPU experiment')
    pairs = packed.pivot(index=['case', 'block'], columns='kernel', values=['time_s', 'energy_j'])
    gains = pairs.xs('direct', level=1, axis=1) / pairs.xs('packed', level=1, axis=1)
    return depths, gains.reset_index()


def make_figure(root: Path):
    depths, gains = read_trials(root)
    cases = [(region, length) for length in (65536, 524288)
             for region in ('collinear', 'insertion')]
    rows = [4, 8, 16, 32]
    labels = [f'{length:,} nt   ·   {"colinéaire" if region == "collinear" else "insertion"}'
              .replace(',', ' ') for region, length in cases]
    use_style()
    figure, (matrix, ratios) = plt.subplots(1, 2, figsize=(8.6, 3.35), width_ratios=[1.35, 1],
                                           layout='constrained')
    records = {}
    extremes = [1.]
    for i, (region, length) in enumerate(cases):
        case = f'{region}-{length}'
        subset = depths[depths.case == case]
        for j, depth in enumerate(rows):
            records[i, j] = sample(subset.loc[subset.rows == depth, 'gcups_w'])
        for key, dy, color, marker, label in (
                ('time_s', -.14, COLORS[0], 'o', 'Temps'),
                ('energy_j', .14, COLORS[1], 'D', 'Énergie')):
            result = sample(gains.loc[gains.case == case, key])
            estimate(ratios, i + dy, result, color, marker=marker, horizontal=True,
                     label=label if not i else None)
            extremes.extend([result['low95'], result['high95']])
    measurement_matrix(matrix, records, labels, list(map(str, rows)), precision='.3f',
                       spread='.2g', color=2)
    ratios.axvline(1, **GUIDE)
    span = np.ptp(extremes)
    ratios.set(xlim=(min(extremes) - .06*span, max(extremes) + .06*span),
               ylim=matrix.get_ylim(), yticks=range(4), yticklabels=['']*4)
    matrix.set_xlabel('Profondeur de tuile [lignes]')
    ratios.set_xlabel('Rapport direct / compacté [1]')
    label_panels((matrix, ratios))
    column_titles((matrix, ratios), ('Rendement [GCUPS/W]', 'Indices de dernière occurrence'))
    legend_above(figure, (ratios,), 2)
    return figure
