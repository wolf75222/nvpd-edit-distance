"""Aggregate the randomized energy trials and compare tile depths within blocks."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'python'))
from bench_stats import sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    data = pd.read_csv(args.directory / 'trials.csv')
    refs = pd.read_csv(ROOT / 'results/roofline/gpu/839191/trials.csv')
    cases = {f'{region}-{size}' for region in ('collinear', 'insertion')
             for size in (65536, 524288)}
    oracle = refs[refs.method.isin(['edlib', 'wfa2']) & refs.case.isin(cases)]
    by_case = oracle.groupby('case')
    if (set(oracle.case) != cases or not by_case.method.nunique().eq(2).all()
            or not by_case.distance.nunique().eq(1).all()):
        raise ValueError('missing score agreement between edlib and WFA2')
    checked = data.merge(by_case.distance.first().rename('reference'), on='case', validate='m:1')
    groups = checked.groupby(['case', 'rows']).block
    if (len(checked) != 80 or len(checked) != len(data) or groups.ngroups != 16
            or not groups.nunique().eq(5).all() or set(checked.block) != set(range(1, 6))
            or set(checked.rows) != {4, 8, 16, 32} or not checked.distance.eq(checked.reference).all()
            or not np.isfinite(checked[['cold_s', 'warm_s', 'seconds', 'joules']]).all().all()
            or not checked[['iterations', 'seconds', 'joules']].gt(0).all().all()):
        raise ValueError('incomplete or inconsistent energy trials')
    checked['time_s'] = checked.seconds / checked.iterations
    checked['energy_j'] = checked.joules / checked.iterations
    checked['power_w'] = checked.joules / checked.seconds
    checked['gcups'] = checked.m * checked.n / (checked.time_s * 1e9)
    checked['gcups_per_w'] = checked.m * checked.n / (checked.energy_j * 1e9)
    base = checked[checked.rows == 8].set_index(['case', 'block'])[['time_s', 'energy_j']]
    checked = checked.join(base, on=['case', 'block'], rsuffix='_8', validate='m:1')
    checked['speedup_vs8'] = checked.time_s_8 / checked.time_s
    checked['energy_gain_vs8'] = checked.energy_j_8 / checked.energy_j
    metrics = ['cold_s', 'warm_s', 'time_s', 'energy_j', 'power_w', 'gcups', 'gcups_per_w',
               'speedup_vs8', 'energy_gain_vs8']
    stats = []
    for (case, rows), group in checked.groupby(['case', 'rows']):
        for metric in metrics:
            estimate = sample(group[metric])
            stats.append(dict(case=case, rows=rows, metric=metric, **estimate))
    pd.DataFrame(stats).to_csv(args.directory / 'summary.csv', index=False, float_format='%.10g')


if __name__ == '__main__':
    main()
