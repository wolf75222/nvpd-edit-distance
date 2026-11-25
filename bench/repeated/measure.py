#!/usr/bin/env python3
"""Run matched, independently restarted blocks with a full same-input warmup."""
import argparse
import csv
import hashlib
import math
import os
from pathlib import Path
import platform
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
FIELDS = ('job platform host lane cpu_ids threads case algorithm block order warmup m n distance seconds '
          'process_peak_rss_kib device_bytes input_sha256 binary_sha256').split()


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def executable(algorithm):
    return "repeat_initial" if algorithm == "nvpd-initial" else "repeat_kernel"


def digest(path):
    d = hashlib.sha256()
    with path.open('rb') as source:
        for part in iter(lambda: source.read(1 << 20), b''):
            d.update(part)
    return d.hexdigest()


def plan(machine, lane):
    if machine=='epyc-scalar':
        yield 'syn-1000000-' + ('0.02','0.05','0.1')[lane % 3], ['seq'], 1
        return
    if machine=='epyc-trace':
        for size in (2000, 8000, 25000):
            yield f'prefix-{size}', ['seq-trace', 'hirschberg-trace'] if size < 25000 else ['hirschberg-trace'], 31
        return
    rates = ('0.02', '0.05', '0.1')
    if machine == 'epyc':
        rate = rates[lane % 3]
        if lane < 3:
            yield 'syn-1000000-' + rate, ['nvpd-initial', 'nvpd-simd'], 5
        else:
            for n in (5000, 20000, 100000):
                yield f'syn-{n}-{rate}', ['seq', 'diagonal', 'nvpd-simd', 'edlib', 'wfa2'], 9 if n < 100000 else 5
            yield 'syn-1000000-' + rate, ['diagonal', 'edlib', 'wfa2'], 5
        return
    for n in (5000, 20000, 100000, 1000000):
        for rate in rates:
            yield f'syn-{n}-{rate}', ['nvpd-cuda', 'edlib', 'wfa2'], 9 if n < 100000 else 5
    for name in ('chr20-assemblies', 'chr21-full-insert64', 'chr21-2m', 'chr21-4m'):
        methods = ['edlib', 'wfa2']
        if name.startswith(('chr21-2m','chr21-4m')):
            methods.append('nvpd-cuda')
        elif name == 'chr21-full-insert64':
            methods.append('banded-trim')
        yield name, methods, 5


def order(methods, block):
    # Cyclic order, reversing on alternate cycles, balances first position.
    k = (block - 1) % len(methods)
    base = methods if (block - 1) // len(methods) % 2 == 0 else methods[::-1]
    return base[k:] + base[:k]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('platform', choices=['epyc','gh200','epyc-trace','epyc-scalar'])
    p.add_argument('build', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--lane', type=int, default=0)
    args = p.parse_args()
    lanes = {'epyc': 6, 'epyc-scalar': 15, 'epyc-trace': 1, 'gh200': 1}
    if not 0 <= args.lane < lanes[args.platform]:
        p.error(f"lane must be between 0 and {lanes[args.platform] - 1} for {args.platform}")
    cases = {row['case']: row for row in read_csv(HERE / 'inputs.csv')}
    common = dict(job=os.environ.get('SLURM_JOB_ID','local'), platform=args.platform,
                  host=platform.node().split('.')[0],
                  lane=args.lane, cpu_ids=';'.join(map(str,sorted(os.sched_getaffinity(0)))) if hasattr(os,'sched_getaffinity') else '',
                  threads=int(os.environ.get('OMP_NUM_THREADS','1')))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    if args.platform=='epyc-scalar':
        cpu_ids=sorted(os.sched_getaffinity(0))
        if len(cpu_ids)!=1 or common['threads']!=1:
            raise ValueError('expected one CPU per scalar process')
        topology=Path(f'/sys/devices/system/cpu/cpu{cpu_ids[0]}/topology')
        place=[str(cpu_ids[0])]+[(topology/k).read_text().strip() for k in ('physical_package_id','core_id')]
        args.out.with_suffix('.affinity').write_text(','.join(place)+'\n')
    hashes = {}
    with args.out.open('x',newline='') as out, tempfile.TemporaryDirectory() as tmp:
        w=csv.DictWriter(out,fieldnames=FIELDS,lineterminator='\n')
        w.writeheader()
        for name, methods, nblocks in plan(args.platform,args.lane):
            case=cases[name]
            path = HERE.parents[1] / case['path']
            if name.startswith('prefix-'):
                source = (HERE.parents[1] / cases['syn-100000-0.02']['path']).read_bytes()
                records=[b''.join(r.splitlines()[1:])[:int(case['m'])] for r in source.split(b'>')[1:]]
                path=Path(tmp)/path.name
                path.write_bytes(b'>A\n'+records[0]+b'\n>B\n'+records[1]+b'\n')
            if digest(path)!=case['input_sha256']:
                raise ValueError(f'input hash mismatch for {name}')
            print(f'lane {args.lane}: {name}, {nblocks} blocks',flush=True)
            blocks = [args.lane//3+1] if args.platform=='epyc-scalar' else range(1,nblocks+1)
            for block in blocks:
                for pos,alg in enumerate(order(methods,block),1):
                    binpath = args.build / executable(alg)
                    if binpath not in hashes:
                        hashes[binpath]=digest(binpath)
                    rows=subprocess.check_output([binpath.resolve(),path,alg,case['distance']],
                                                 text=True,input='run\nrun\n').splitlines()
                    if len(rows)!=3 or rows[0]!='ready':
                        raise ValueError(f'unexpected result for {name} {alg}: {rows}')
                    for warm, line in zip((1,0),rows[1:]):
                        m,n,dist,sec,rss,dev=line.split(',')
                        seconds = float(sec)
                        if ([m,n,dist]!=[case[k] for k in ('m','n','distance')] or
                                not math.isfinite(seconds) or seconds <= 0 or int(rss) < 0 or int(dev) < 0):
                            raise ValueError(f'result mismatch for {name} {alg}')
                        w.writerow(dict(**common,case=name,algorithm=alg,block=block,order=pos,
                                        warmup=warm,m=m,n=n,distance=dist,seconds=sec,process_peak_rss_kib=rss,
                                        device_bytes=dev,input_sha256=case['input_sha256'],binary_sha256=hashes[binpath]))
                    out.flush()


if __name__ == '__main__':
    main()
