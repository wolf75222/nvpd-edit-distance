#!/usr/bin/env python3
"""Validate every warmup and timed block, then publish per-input Student intervals."""
import argparse
import csv
from collections import defaultdict
from itertools import combinations
from math import isfinite, sqrt
from pathlib import Path
from statistics import mean, median, stdev

from measure import HERE, FIELDS, order, plan

# Two-sided 95% t critical values for the experiment's fixed sample sizes.
T95 = {5: 2.7764451051977987, 8: 2.364624251592785, 9: 2.306004135204166, 31: 2.0422724563012373}


def write(path, rows):
    if not rows:
        return
    tmp = path.with_suffix('.tmp')
    with tmp.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=rows[0],lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def stats(data):
    n = len(data)
    mu, sd = mean(data), stdev(data)
    half = T95[n]*sd/sqrt(n)
    return dict(blocks=n,mean=mu,sd=sd,ci95_low=mu-half,ci95_high=mu+half,
                median=median(data),ri95=half/mu,within_5pct=int(half/mu<=0.05))


def memory_stats(peaks):
    # One terminal high-water mark per independently started process.
    valid=[v for v in peaks if v>0]
    result=dict(rss_observed_processes=len(valid),max_process_peak_rss_kib=max(valid,default=''))
    for key in ('mean','sd','ci95_low','ci95_high'):
        result[key+'_process_peak_rss_kib']=''
    if len(valid)==len(peaks):
        estimate=stats(valid)
        for key in ('mean','sd','ci95_low','ci95_high'):
            result[key+'_process_peak_rss_kib']=estimate[key]
    return result


def collect_scalar(raw, output, cases):
    trials, placement, summary = [], [], []
    for lane in range(15):
        case = next(plan('epyc-scalar',lane))[0]
        ref = cases[case]
        rows = list(csv.DictReader((raw/f'lane-{lane}.csv').open()))
        cpu,socket,core=(raw/f'lane-{lane}.affinity').read_text().strip().split(',')
        if len(rows)!=2 or [r['warmup'] for r in rows]!=['1','0']:
            raise ValueError(f'incomplete scalar lane {lane}')
        for r in rows:
            if not (r['case']==case and r['algorithm']=='seq' and int(r['lane'])==lane and
                    r['platform']=='epyc-scalar' and int(r['block'])==lane//3+1 and
                    r['threads']=='1' and r['cpu_ids']==cpu and r['order']=='1' and
                    all(r[k]==ref[k] for k in ('m','n','distance','input_sha256')) and
                    isfinite(float(r['seconds'])) and float(r['seconds'])>0):
                raise ValueError(f'inconsistent scalar lane {lane}')
        trials.extend(rows)
        placement.append(dict(lane=lane,case=case,block=lane//3+1,cpu=cpu,socket=socket,physical_core=core))
    if len({(r['socket'],r['physical_core']) for r in placement})!=15:
        raise ValueError('scalar tasks share a physical core')
    if any(len({r[k] for r in trials})!=1 for k in ('host','job','binary_sha256')):
        raise ValueError('scalar workers differ in placement, allocation or binary')
    for case in sorted({r['case'] for r in trials}):
        group=[r for r in trials if r['case']==case and r['warmup']=='0']
        if {int(r['block']) for r in group}!=set(range(1,6)):
            raise ValueError(f'incomplete scalar processes for {case}')
        summary.append(dict(platform='epyc-scalar',case=case,input_sha256=cases[case]['input_sha256'],
                            job=trials[0]['job'],experimental_unit='fresh_process_on_distinct_core_after_one_warmup',
                            algorithm='seq',cpu_workers=1,**stats([float(r['seconds']) for r in group]),
                            **memory_stats([int(r['process_peak_rss_kib']) for r in group])))
    output.mkdir(parents=True,exist_ok=True)
    write(output/'trials.csv',trials)
    write(output/'summary.csv',summary)
    write(output/'placement.csv',placement)
    print('EPYC scalar: 15 timed calls, 15 warmups on 15 distinct physical cores.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('platform',choices=['epyc','gh200','epyc-trace','epyc-scalar'])
    p.add_argument('raw',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    cases={r['case']:r for r in csv.DictReader((HERE/'inputs.csv').open())}
    if args.platform=='epyc-scalar':
        collect_scalar(args.raw,args.output,cases)
        return
    trials=[]
    summary=[]
    ratios=[]
    placement=set()
    for lane in range(6 if args.platform=='epyc' else 1):
        rows=list(csv.DictReader((args.raw/f'lane-{lane}.csv').open()))
        if list(rows[0])!=FIELDS:
            raise ValueError('unexpected raw columns')
        masks={r['cpu_ids'] for r in rows}
        if len(masks)!=1:
            raise ValueError('placement changed during a lane')
        used=set(next(iter(masks)).split(';'))
        expected_cores={'epyc':8,'epyc-trace':1,'gh200':4}.get(args.platform)
        if expected_cores is not None and len(used)!=expected_cores:
            raise ValueError(f'unexpected CPU affinity mask {masks}')
        if placement & used:
            raise ValueError('overlapping CPU lanes')
        placement.update(used)
        lookup={(r['case'],r['algorithm'],int(r['block']),int(r['warmup'])):r for r in rows}
        if len(lookup)!=len(rows):
            raise ValueError('duplicate result')
        for name,methods,count in plan(args.platform,lane):
            ref=cases[name]
            times=defaultdict(list)
            peaks=defaultdict(list)
            for block in range(1,count+1):
                for pos,alg in enumerate(order(methods,block),1):
                    for warm in (1,0):
                        r=lookup.pop((name,alg,block,warm))
                        matches=all(r[k]==ref[k] for k in ('m','n','distance','input_sha256'))
                        t=float(r['seconds'])
                        if not (matches and r['platform']==args.platform and int(r['lane'])==lane
                                and int(r['order'])==pos and int(r['threads'])==(8 if args.platform=='epyc' else 1)
                                and int(r['process_peak_rss_kib'])>=0 and int(r['device_bytes'])>=0
                                and isfinite(t) and t>0):
                            raise ValueError(f'invalid observation {name} {alg} {block}')
                        trials.append(r)
                        if not warm:
                            times[alg].append(t)
                            peaks[alg].append(int(r['process_peak_rss_kib']))
            common=dict(platform=args.platform,case=name,input_sha256=ref['input_sha256'],lane=lane,
                        job=rows[0]['job'],experimental_unit='fresh_process_after_one_warmup')
            for alg in methods:
                record=stats(times[alg])
                workers=8 if args.platform=='epyc' and alg in ('nvpd-initial','nvpd-simd','diagonal') else 1
                summary.append(dict(**common,algorithm=alg,cpu_workers=workers,**record,
                                    **memory_stats(peaks[alg])))
            for a,b in combinations(methods,2):
                values=[u/v for u,v in zip(times[a],times[b])]
                ratios.append(dict(**common,numerator=a,denominator=b,**stats(values)))
        if lookup:
            raise ValueError('unexpected observations')
    if len({r['job'] for r in trials})!=1 or len({r['host'] for r in trials})!=1:
        raise ValueError('multiple allocations or hosts')
    args.output.mkdir(parents=True,exist_ok=True)
    write(args.output/'trials.csv',trials)
    write(args.output/'summary.csv',summary)
    write(args.output/'paired.csv',ratios)
    print(f'{args.platform}: {len(trials)//2} measured calls + same number of warmups, '
          f'{len(summary)} cells, {sum(r["within_5pct"] for r in summary)} within 5% relative half-width')


if __name__=='__main__':
    main()
