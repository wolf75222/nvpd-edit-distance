#!/usr/bin/env python3
"""Match the QuickEd EPYC protocol on the two banded algorithms and library controls."""
import argparse
import csv
from itertools import combinations
import os
from pathlib import Path
import random
from statistics import mean
import subprocess

from collect import memory_stats, stats, write
from measure import digest

METHODS=['banded','banded-trim','edlib','wfa2']
CASES=['20k-0.02','20k-0.05','20k-0.1','chr21-50k','chr21-center128']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary',type=Path)
    p.add_argument('inputs',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    cases={r['case']:r for r in csv.DictReader((args.inputs/'cases.csv').open())}
    affinity=sorted(os.sched_getaffinity(0))
    if len(affinity)!=1 or os.environ['OMP_NUM_THREADS']!='1':
        raise ValueError('one pinned core required')
    rng=random.Random(2026)
    trials,summary,ratios=[],[],[]
    args.output.mkdir(parents=True,exist_ok=True)
    h=digest(args.binary)
    for name in CASES:
        case=cases[name]
        path=args.inputs/(name+'.fa')
        if digest(path)!=case['sha256']:
            raise ValueError(f'hash mismatch {name}')
        common=dict(platform='epyc-bands',case=name,input_sha256=case['sha256'])
        averages={alg:[] for alg in METHODS}
        peaks={alg:[] for alg in METHODS}
        for block in range(1,9):
            methods=METHODS.copy()
            rng.shuffle(methods)
            for pos,alg in enumerate(methods,1):
                lines=subprocess.check_output([args.binary,path,alg,case['exact']],text=True,input='run\n'*11).splitlines()
                if lines[0]!='ready' or len(lines)!=12:
                    raise ValueError(f'missing measurements {name} {alg}')
                measured=[]
                for call,line in enumerate(lines[1:]):
                    m,n,score,seconds,rss,dev=line.split(',')
                    if [m,n,score]!=[case[k] for k in ('m','n','exact')] or not float(seconds)>0:
                        raise ValueError(f'invalid result {name} {alg}')
                    trials.append(dict(**common,algorithm=alg,block=block,order=pos,call=call,
                                       warmup=int(call==0),m=m,n=n,distance=score,seconds=seconds,
                                       process_peak_rss_kib=rss,device_bytes=dev,job=os.environ['SLURM_JOB_ID'],
                                       cpu_id=affinity[0],binary_sha256=h))
                    if call:
                        measured.append(float(seconds))
                averages[alg].append(mean(measured))
                peaks[alg].append(int(rss))
        for alg in METHODS:
            summary.append(dict(**common,algorithm=alg,cpu_workers=1,**stats(averages[alg]),**memory_stats(peaks[alg])))
        for a,b in combinations(METHODS,2):
            ratios.append(dict(**common,numerator=a,denominator=b,
                               **stats([u/v for u,v in zip(averages[a],averages[b])])))
        print(f'{name}: 8 blocks x 4 methods x 10 checked measurements',flush=True)
    for row in summary + ratios:
        row.update(job=os.environ['SLURM_JOB_ID'],
                   experimental_unit='mean_of_ten_calls_in_one_fresh_warmed_process')
    write(args.output/'trials.csv',trials)
    write(args.output/'summary.csv',summary)
    write(args.output/'paired.csv',ratios)


if __name__=='__main__':
    main()
