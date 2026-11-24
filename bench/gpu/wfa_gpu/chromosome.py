#!/usr/bin/env python3
"""Measure the GPU/CPU handoff threshold on homologous human chromosome windows."""
import csv
import hashlib
import random
import sys
from pathlib import Path

from measure import BLOCKS, ROOT, estimates, run


def prepare(fasta, destination):
    records = {}
    for line in fasta.read_bytes().splitlines():
        if line.startswith(b'>'):
            name = line[1:].decode()
            records[name] = bytearray()
        else:
            records[name].extend(line)
    base = ROOT / 'results/mpi/human'
    with (base / 'scores.tsv').open() as source:
        scores = {(r['id_a'],r['id_b']): int(r['edlib']) for r in csv.DictReader(source, delimiter='\t')}
    expected = []
    with (base / 'windows.tsv').open() as source, destination.open('wb') as output:
        for row in csv.DictReader(source, delimiter='\t'):
            if not row['id'].startswith('locus'):
                continue
            keys = [row['id'] + suffix for suffix in ('_A','_B')]
            for key, column, marker in zip(keys, ('sha256_a','sha256_b'), (b'>',b'<')):
                sequence = records[key]
                if hashlib.sha256(sequence).hexdigest() != row[column]:
                    raise ValueError(f'{key}: unexpected sequence digest')
                output.write(marker + sequence + b'\n')
            expected.append(scores[tuple(keys)])
    if len(expected) != 768:
        raise ValueError('expected 768 homologous pairs')
    return expected


def measure(work, out, fasta):
    path = work / 'human.seq'
    expected = prepare(fasta, path)
    count = len(expected)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    configs = [('edlib',0), ('wfa2',0), ('nvpd-cuda',0)] + [('wfa-gpu',k) for k in (64,256,1024)]
    times = {c: [] for c in configs}
    rng = random.Random(20261124)
    with (out / 'trials.csv').open('w') as file:
        csvout = csv.writer(file)
        csvout.writerow('case,batch_sha256,pairs,cap,block,position,method,seconds,cpu_fallback,wrong_scores,min_expected,max_expected'.split(','))
        for block in range(BLOCKS + 1):
            order = configs.copy()
            rng.shuffle(order)
            for position, (method, cap) in enumerate(order):
                elapsed, values, fallbacks = run(method, work, path, count, cap)
                wrong = sum(a != b for a,b in zip(values, expected))
                csvout.writerow(('chr21',digest,count,cap,block,position,method,elapsed,fallbacks,wrong,min(expected),max(expected)))
                file.flush()
                if wrong:
                    raise RuntimeError(f'{wrong} wrong global scores for {method} at cap {cap}')
                if block:
                    times[method,cap].append(elapsed)
            print('chr21 block', block, 'completed', flush=True)
    with (out / 'summary.csv').open('w') as file:
        wr = csv.writer(file)
        wr.writerow('case,pairs,method,cap,blocks,mean_s,sd_s,ci95_low_s,ci95_high_s,ratio_over_wfa2,ratio_sd,ratio_lo,ratio_hi'.split(','))
        for (method,cap), values in times.items():
            ratios = [v/w for v,w in zip(values, times['wfa2',0])]
            wr.writerow(('chr21',count,method,cap,BLOCKS,*estimates(values),*estimates(ratios)))


if __name__ == '__main__':
    measure(*map(Path, sys.argv[1:]))
