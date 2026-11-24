#!/usr/bin/env python3
"""Build a memory series without changing any GH200 timing or combining distinct metrics."""
import argparse
import csv
from pathlib import Path
import sys
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from collect import stats, write


def group(path):
    grouped=defaultdict(list)
    for row in csv.DictReader(path.open()):
        grouped[row['case'],row['algorithm']].append(row)
    return grouped


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('gh200',type=Path)
    p.add_argument('probe',type=Path)
    args=p.parse_args()
    primary=group(args.gh200/'trials.csv')
    new=group(args.probe/'trials.csv')
    trials,summary=[],[]
    for (case,alg),calls in primary.items():
        calls=[r for r in calls if r['warmup']=='0']
        chosen=calls
        metric='getrusage_after_timed_call'
        key='process_peak_rss_kib'
        count=len(calls)
        if not all(int(r[key])>0 for r in calls):
            chosen=new.pop((case,alg))
            by_block={x['block']:x for x in calls}
            if len(chosen)!=count or {r['block'] for r in chosen}!=set(by_block):
                raise ValueError(f'incomplete memory series for {case} {alg}')
            for r in chosen:
                gold=by_block[r['block']]
                if any(r[k]!=gold[k] for k in ('input_sha256','distance','m','n')) or r['worker_cpu']==r['observer_cpu']:
                    raise ValueError('input, score or affinity mismatch')
                if int(r['samples'])<5 or int(r['sampled_smaps_rss_kib'])<=0:
                    raise ValueError('missing RSS samples or barriers')
            metric='sampled_smaps_rollup_max'
            key='sampled_smaps_rss_kib'
            for native,label in [('wait4_maxrss_kib','final_parent_wait4'),('status_vmhwm_kib','barrier_status_VmHWM')]:
                if all(int(r[native])>0 and int(r[native])>=int(r['sampled_smaps_rss_kib']) for r in chosen):
                    key,metric=native,label
                    break
        peaks=[int(r[key]) for r in chosen]
        info=dict(case=case,algorithm=alg,input_sha256=calls[0]['input_sha256'],timing_job=calls[0]['job'],
                  memory_job=chosen[0]['job'],rss_metric=metric)
        for r,value in zip(chosen,peaks):
            trials.append(dict(**info,block=r['block'],rss_kib=value))
        st=stats(peaks)
        summary.append(dict(**info,processes=count,mean_rss_kib=st['mean'],sd_rss_kib=st['sd'],
                            ci95_low_rss_kib=st['ci95_low'],ci95_high_rss_kib=st['ci95_high'],
                            max_rss_kib=max(peaks),experimental_unit='terminal_or_sampled_maximum_for_one_warmed_process'))
    if new or len(summary)!=47:
        raise ValueError('memory coverage differs from the 47 timing configurations')
    write(args.gh200/'memory_trials.csv',trials)
    write(args.gh200/'memory.csv',summary)
    print(f'{len(trials)} per-process memory observations for {len(summary)} configurations.')


if __name__=='__main__':
    main()
