#!/usr/bin/env python3
"""Score and time a literal-IUPAC cohort from FASTA to pairwise TSV."""
import csv
import hashlib
import itertools
import random
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from pathlib import Path

from measure import BLOCKS, ROOT, estimates

WORKERS = 4
CAP = 256
METHODS = ('wfa2', 'edlib', 'hybrid')


def records(fasta):
    seqs = {}
    for line in fasta.read_bytes().splitlines():
        if line.startswith(b'>'):
            key = line[1:].split()[0].decode('ascii')
            if key in seqs:
                raise ValueError(f'duplicate accession {key}')
            seqs[key] = bytearray()
        elif line:
            seqs[key].extend(line.strip().upper())
    return seqs


def score(work, fasta, method):
    seqs = records(fasta)
    eligible = {key for key,s in seqs.items() if len(s) < 32768 and set(s).issubset(b'ACGT')}
    pairs = list(itertools.combinations(seqs, 2))
    shards = [[] for _ in range(WORKERS + 1)]
    with ExitStack() as stack:
        streams = [stack.enter_context((work / f'shard{i}.seq').open('wb')) for i in range(len(shards))]
        for idx, (a,b) in enumerate(pairs):
            use_gpu = method == 'hybrid' and a in eligible and b in eligible
            shard = WORKERS if use_gpu else idx % WORKERS
            shards[shard].append(idx)
            streams[shard].write(b'>' + seqs[a] + b'\n<' + seqs[b] + b'\n')

    def compute(shard):
        path, output = work / f'shard{shard}.seq', work / f'shard{shard}.tsv'
        gpu = shard == WORKERS
        if gpu:
            cmd = [work / 'upstream/bin/wfa.affine.gpu', '-g', '1,0,1', '-e', str(CAP), '-b',
                   str(min(512, len(shards[shard]))), '-i', path, '-o', output]
        else:
            backend = 'wfa2' if method == 'hybrid' else method
            cmd = [work / 'baselines/wfagpu_baselines', backend, path, output]
        start = time.perf_counter()
        result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, check=True)
        elapsed = time.perf_counter() - start
        sign = -1 if gpu else 1
        scores = [sign * int(line.split()[0]) for line in output.read_text().splitlines()]
        if len(scores) != len(shards[shard]):
            raise ValueError(f'incomplete output from {method} shard {shard}')
        fallback = sum(map(int, re.findall(r'\) (\d+)/\d+ alignemnts could not', result.stderr)))
        return shard, elapsed, fallback, scores

    scores = [None] * len(pairs)
    longest_cpu, gpu_elapsed, fallback = 0., 0., 0
    with ThreadPoolExecutor(max_workers=len(shards)) as pool:
        for shard, elapsed, fallback_count, values in pool.map(compute, [i for i,s in enumerate(shards) if s]):
            if shard == WORKERS:
                gpu_elapsed, fallback = elapsed, fallback_count
            else:
                longest_cpu = max(longest_cpu, elapsed)
            for index,value in zip(shards[shard], values):
                scores[index] = value
    with (work / 'cohort_scores.tsv').open('w') as dest:
        wr = csv.writer(dest, delimiter='\t', lineterminator='\n')
        wr.writerow(('id_a','id_b','distance'))
        for (a,b), value in zip(pairs, scores):
            wr.writerow((a,b,value))
    print(len(shards[-1]), fallback, longest_cpu, gpu_elapsed)


def measure(work, out, fasta):
    inputs = records(fasta)
    with (ROOT / 'results/ncbi/sars_content.tsv').open() as source:
        manifest = {r['id']: r['sha256'] for r in csv.DictReader(source, delimiter='\t')}
    if len(inputs) != 64 or {k: hashlib.sha256(s).hexdigest() for k,s in inputs.items()} != manifest:
        raise ValueError('the 64 records must match the SARS-CoV-2 project manifest')
    with (ROOT / 'results/phylogeny/sars-edit/matrix.tsv').open() as source:
        matrix = {r['id']: r for r in csv.DictReader(source, delimiter='\t')}
    reference = [(a,b,int(matrix[a][b])) for a,b in itertools.combinations(inputs,2)]
    times = {m: [] for m in METHODS}
    rng = random.Random(20261125)
    with (out / 'trials.csv').open('w') as file:
        wr = csv.writer(file)
        wr.writerow('block,position,method,cpu_workers,seconds,gpu_submitted,gpu_fallback,cpu_scored,gpu_scored,max_cpu_shard_s,gpu_shard_s,pairs,wrong_scores'.split(','))
        for block in range(BLOCKS + 1):
            order = list(METHODS)
            rng.shuffle(order)
            for position, method in enumerate(order):
                start = time.perf_counter()
                child = subprocess.run([sys.executable,__file__,'score',work,fasta,method], capture_output=True, text=True)
                elapsed = time.perf_counter() - start
                if child.returncode:
                    raise RuntimeError(child.stderr)
                submitted, fallback, cpu_time, gpu_time = child.stdout.split()
                active = int(submitted) - int(fallback)
                with (work / 'cohort_scores.tsv').open() as source:
                    current = [(r['id_a'],r['id_b'],int(r['distance'])) for r in csv.DictReader(source, delimiter='\t')]
                if current != reference:
                    raise ValueError(f'{method} did not reproduce the literal matrix')
                wr.writerow((block,position,method,WORKERS,elapsed,submitted,fallback,len(reference)-active,active,cpu_time,gpu_time,len(reference),0))
                file.flush()
                if block:
                    times[method].append(elapsed)
            print('SARS-CoV-2 block', block, 'completed', flush=True)
    with (out / 'summary.csv').open('w') as file:
        wr = csv.writer(file)
        wr.writerow('method,blocks,pairs,mean_s,sd_s,ci95_low_s,ci95_high_s,ratio_over_wfa2,ratio_sd,ratio_lo,ratio_hi'.split(','))
        for method, values in times.items():
            wr.writerow((method,BLOCKS,len(reference),*estimates(values),*estimates([t/w for t,w in zip(values,times['wfa2'])])))


if __name__ == '__main__':
    if sys.argv[1] == 'score':
        score(Path(sys.argv[2]),Path(sys.argv[3]),sys.argv[4])
    else:
        measure(*map(Path, sys.argv[1:]))
