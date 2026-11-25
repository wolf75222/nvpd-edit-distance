#!/usr/bin/env python3
"""Summarize paired ratios across independently generated sequence pairs."""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'python'))
import ct_plot


def analyze(path):
    data = pd.read_csv(path)
    if len(data) != 1080 or data.duplicated(['path', 'block', 'algorithm', 'warmup']).any():
        raise ValueError('Expected one warm-up and one timing per pair, block and method')
    if data.groupby('path')[['sha256', 'distance']].nunique().ne(1).any().any():
        raise ValueError('An input changed or produced different scores')
    if not np.isfinite(data.seconds).all() or (data.seconds <= 0).any():
        raise ValueError('Invalid execution time')
    measured = data[data.warmup == 0]
    keys = ['length', 'rate', 'rep', 'block']
    times = measured.pivot(index=keys, columns='algorithm', values='seconds')
    if len(times) != 60 * 3 or times.isna().any().any():
        raise ValueError('Expected three paired blocks for each of 60 inputs')
    ratios = pd.DataFrame({name: times['nvpd-cuda'] / times[name] for name in ('edlib', 'wfa2')})
    pairs = ratios.groupby(level=['length', 'rate', 'rep']).mean().stack().rename('ratio').reset_index()
    pairs = pairs.rename(columns={'level_3': 'reference'})
    groups = pairs.groupby(['length', 'rate', 'reference']).ratio
    result = groups.agg(['mean', 'std', 'count']).reset_index()
    if not result['count'].eq(5).all():
        raise ValueError('Expected five independent input pairs per configuration')
    result['half95'] = t.ppf(.975, 4) * result['std'] / np.sqrt(5)
    return pairs, result


def plot(pairs, summary):
    ct_plot.use_style()
    fig, axes = plt.subplots(1, 3, figsize=(8.0, 3.0), sharex=True, sharey=True)
    lengths = [5000, 20000, 100000, 1000000]
    labels = ['5', '20', '100', '1 000']
    for ax, rate in zip(axes, [.02, .05, .1]):
        for j, reference in enumerate(('edlib', 'wfa2')):
            rows = summary[(summary.rate == rate) & (summary.reference == reference)].set_index('length').loc[lengths]
            positions = np.arange(4) + (j - .5) * .18
            ax.errorbar(positions, rows['mean'], yerr=rows['half95'], fmt='o', markersize=4,
                        capsize=2.5, color=ct_plot.COLORS[j], label='Edlib' if j == 0 else 'WFA2')
            ax.errorbar(positions, rows['mean'], yerr=rows['std'], fmt='none',
                        color=ct_plot.COLORS[j], elinewidth=3, alpha=.55)
        ax.axhline(1, color=ct_plot.INK, linewidth=.8, linestyle='--')
        ax.set_yscale('log')
        ax.set_xticks(range(4), labels)
        ax.set_xlabel('Longueur [knt]')
    axes[0].set_ylabel('Temps NvPD / temps de référence')
    ct_plot.column_titles(axes, ['2 %', '5 %', '10 %'])
    ct_plot.label_panels(axes)
    ct_plot.legend_above(fig, axes, ncol=2)
    return fig


def main():
    p = argparse.ArgumentParser()
    p.add_argument('trials', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--figure', type=Path, required=True)
    args = p.parse_args()
    pairs, summary = analyze(args.trials)
    summary.to_csv(args.output, index=False)
    figure = plot(pairs, summary)
    figure.savefig(args.figure)
    plt.close(figure)


if __name__ == '__main__':
    main()
