#!/usr/bin/env python3
"""Measure per-process memory independently of the timing campaign on Grace."""
import argparse
import csv
import os
from pathlib import Path
import sys
import select
import subprocess
from time import monotonic

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from measure import HERE, digest

POLL_S = .001


def proc_kib(pid, file, key):
    with Path(f'/proc/{pid}/{file}').open() as stream:
        for line in stream:
            if line.startswith(key+':'):
                value,unit=line.split()[1:]
                if unit!='kB':
                    raise ValueError(f'unexpected proc unit {unit}')
                return int(value)
    raise ValueError(f'{key} missing from {file}')


def one(binary, path, case, alg, cpu):
    with subprocess.Popen(['taskset','-c',str(cpu),binary,path,alg,case['distance']],
                          stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1) as proc:
        if proc.stdout.readline().strip()!='ready':
            raise RuntimeError('worker did not reach FASTA-loaded barrier')
        vm_peak=child_peak=sampled=interior=0
        samples=0
        max_gap=0.
        last=None

        def sample():
            nonlocal sampled,samples,max_gap,last
            value=proc_kib(proc.pid,'smaps_rollup','Rss')
            now=monotonic()
            if last is not None:
                max_gap=max(max_gap,now-last)
            last=now
            sampled=max(value,sampled)
            samples+=1
            return value

        loaded=sample()
        for _ in range(2):
            vm_peak=max(vm_peak,proc_kib(proc.pid,'status','VmHWM'))
            proc.stdin.write('run\n')
            proc.stdin.flush()
            # Sample immediately, then at nominal 1 ms intervals until the worker
            # signals completion. Very short calls can finish within one interval.
            sample()
            while not select.select([proc.stdout],[],[],POLL_S)[0]:
                sample()
                interior+=1
            m,n,score,_,rss,_=proc.stdout.readline().strip().split(',')
            if [m,n,score]!=[case[k] for k in ('m','n','distance')]:
                raise RuntimeError(f'incorrect score or dimensions: {case["case"]} {alg}')
            child_peak=max(child_peak,int(rss))
            vm_peak=max(vm_peak,proc_kib(proc.pid,'status','VmHWM'))
            sample()
        proc.stdin.close()
        _,status,usage=os.wait4(proc.pid,0)
        proc.returncode=os.waitstatus_to_exitcode(status)
        if proc.returncode or not sampled:
            raise RuntimeError('memory worker failed')
        return dict(loaded_smaps_rss_kib=loaded,wait4_maxrss_kib=usage.ru_maxrss,
                    status_vmhwm_kib=vm_peak,child_maxrss_kib=child_peak,
                    sampled_smaps_rss_kib=sampled,samples=samples,interior_polls=interior,
                    target_poll_ms=POLL_S*1000,max_observed_sample_gap_ms=1000*max_gap)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary',type=Path);p.add_argument('timing_summary',type=Path);p.add_argument('out',type=Path)
    args=p.parse_args()
    cpus=sorted(os.sched_getaffinity(0))
    if len(cpus)<2:
        raise ValueError('observer and child need separate CPUs')
    os.sched_setaffinity(0,{cpus[1]})
    cases={r['case']:r for r in csv.DictReader((HERE/'inputs.csv').open())}
    missing=[r for r in csv.DictReader(args.timing_summary.open()) if not r['mean_process_peak_rss_kib']]
    if len(missing)!=28 or {r['algorithm'] for r in missing}!={'edlib','wfa2'}:
        raise ValueError('unexpected set of missing primary memory observations')
    binary=args.binary.resolve()
    sha=digest(binary)
    args.out.mkdir(parents=True,exist_ok=True)

    def collect(target, experiment):
        with target.open('x',newline='') as stream:
            writer=None
            for name,methods,count in experiment:
                case=cases[name];path=Path(case['path'])
                if digest(path)!=case['input_sha256']:
                    raise ValueError(f'changed FASTA: {name}')
                for block in range(1,count+1):
                    for pos,alg in enumerate(methods if block%2 else methods[::-1],1):
                        metrics=one(binary,path,case,alg,cpus[0])
                        row=dict(job=os.environ['SLURM_JOB_ID'],case=name,algorithm=alg,block=block,
                                 order=pos,input_sha256=case['input_sha256'],binary_sha256=sha,
                                 distance=case['distance'],m=case['m'],n=case['n'],
                                 worker_cpu=cpus[0],observer_cpu=cpus[1],**metrics)
                        if writer is None:
                            writer=csv.DictWriter(stream,fieldnames=row,lineterminator='\n');writer.writeheader()
                        writer.writerow(row);stream.flush()
                print(f'{name}, {count} block(s)',flush=True)
    collect(args.out/'calibration.csv',[('syn-5000-0.02',['edlib','wfa2'],1)])
    for x in csv.DictReader((args.out/'calibration.csv').open()):
        print('calibration',x['algorithm'],*[f'{k}={x[k]}' for k in ('wait4_maxrss_kib','status_vmhwm_kib','sampled_smaps_rss_kib')],flush=True)
    by_case={}
    for r in missing:
        by_case.setdefault(r['case'],[]).append(r)
    collect(args.out/'trials.csv',[(name,[r['algorithm'] for r in rows],int(rows[0]['blocks']))
                                    for name,rows in by_case.items()])


if __name__=='__main__':
    main()
