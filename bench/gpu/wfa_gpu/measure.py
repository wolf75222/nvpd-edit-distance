#!/usr/bin/env python3
"""Wall time and score agreement for fresh processes on matched sequence batches."""
import csv
import hashlib
import itertools
import random
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts/datasets'))
from synthetic import make_pair, pair_seed

METHODS = ('edlib', 'wfa2', 'nvpd-cuda', 'wfa-gpu')
BLOCKS = 5
T95 = 2.7764451051977987  # two-sided t(4)


def batches():
    for length, rate in itertools.product((5000, 20000), ('0.02', '0.1')):
        edits = length * int(float(rate) * 100) // 100
        pairs = [make_pair(length, edits, pair_seed(2026, length, rate, i))[:2]
                 for i in range(1, 129)]
        for count in (5, 128):
            yield f'n{length}_r{rate}_p{count}', pairs[:count], edits + 8, None
    sequences = {}
    for line in (ROOT / 'data/phix174/cohort.fa').read_bytes().splitlines():
        if line.startswith(b'>'):
            accession = line[1:].decode('ascii')
            sequences[accession] = bytearray()
        else:
            sequences[accession].extend(line)
    with (ROOT / 'results/ncbi/phix174.tsv').open() as source:
        manifest = list(csv.DictReader(source, delimiter='\t'))
    if len(sequences) != len(manifest) or len(manifest) != 64:
        raise ValueError('expected 64 phiX174 sequences')
    for row in manifest:
        if hashlib.sha256(sequences[row['id']]).hexdigest() != row['sha256']:
            raise ValueError(f"unexpected sequence for {row['id']}")
    with (ROOT / 'results/romeo/cohort-837271/matrix.tsv').open() as source:
        matrix = {row['id']: row for row in csv.DictReader(source, delimiter='\t')}
    keys = list(itertools.combinations((row['id'] for row in manifest), 2))
    yield 'phix174_all', [(sequences[a], sequences[b]) for a,b in keys], 256, [int(matrix[a][b]) for a,b in keys]


def run(method, work, path, count, threshold):
    scores = work / 'scores.tsv'
    scores.unlink(missing_ok=True)
    if method == 'wfa-gpu':
        cmd = [work / 'upstream/bin/wfa.affine.gpu', '-i', path, '-b', str(min(512, count)),
               '-g', '1,0,1', '-e', str(threshold), '-o', scores]
    else:
        cmd = [work / 'baselines/wfagpu_baselines', method, path, scores]
    start = time.perf_counter()
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    elapsed = time.perf_counter() - start
    if result.returncode or not scores.is_file():
        raise RuntimeError(f'{method} failed: {result.stderr[-2000:]}')
    sign = -1 if method == 'wfa-gpu' else 1
    values = [sign * int(row.split()[0]) for row in scores.read_text().splitlines()]
    if len(values) != count:
        raise RuntimeError(f'{method} returned {len(values)} scores for {count} pairs')
    fallbacks = sum(map(int, re.findall(r'\) (\d+)/\d+ alignemnts could not', result.stderr)))
    return elapsed, values, fallbacks


def estimates(values):
    mean, sd = statistics.mean(values), statistics.stdev(values)
    delta = T95 * sd / len(values) ** 0.5
    return mean, sd, mean - delta, mean + delta


def measure(work, out):
    rng = random.Random(20261124)
    with (out / 'trials.csv').open('w') as trials, (out / 'summary.csv').open('w') as summary:
        tw, sw = csv.writer(trials), csv.writer(summary)
        tw.writerow('case,batch_sha256,pairs,cap,block,position,method,seconds,cpu_fallback,wrong_scores,min_expected,max_expected'.split(','))
        sw.writerow('case,pairs,method,blocks,mean_s,sd_s,ci95_low_s,ci95_high_s,ratio_gpu_over_method,ratio_sd,ratio_lo,ratio_hi'.split(','))
        for name, pairs, threshold, expected in batches():
            path = work / 'pairs.seq'
            with path.open('wb') as target:
                for a,b in pairs:
                    target.write(b'>' + a + b'\n<' + b + b'\n')
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            count = len(pairs)
            if expected is None:
                _, expected, _ = run('edlib', work, path, count, threshold)
            if name == 'n5000_r0.02_p5':
                _, scores, fallbacks = run('wfa-gpu', work, path, count, 32)
                if fallbacks != count or scores != expected:
                    raise RuntimeError(f'fallback control failed: {fallbacks}, {scores}, expected {expected}')
                (out / 'fallback-control.txt').write_text(f'{name} cap=32 cpu={fallbacks}/{count} scores={scores}\n')
            times = {method: [] for method in METHODS}
            for block in range(BLOCKS + 1):
                order = list(METHODS)
                rng.shuffle(order)
                for position, method in enumerate(order):
                    elapsed, values, fallbacks = run(method, work, path, count, threshold)
                    wrong = sum(a != b for a,b in zip(values, expected))
                    if wrong:
                        raise ValueError(f'{name} {method}: {wrong} wrong scores; examples {[(i,v,e) for i,(v,e) in enumerate(zip(values,expected)) if v != e][:8]}')
                    tw.writerow((name,digest,count,threshold,block,position,method,elapsed,fallbacks,wrong,min(expected),max(expected)))
                    trials.flush()
                    if block:
                        times[method].append(elapsed)
            for method in METHODS:
                ratios = [gpu / base for gpu,base in zip(times['wfa-gpu'], times[method])]
                sw.writerow((name,count,method,BLOCKS,*estimates(times[method]),*estimates(ratios)))
            summary.flush()
            print(name, 'completed', flush=True)


if __name__ == '__main__':
    measure(*map(Path, sys.argv[1:]))
