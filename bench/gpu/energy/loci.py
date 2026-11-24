"""Compare two CHM13 candidate loci against one GRCh38 segment through the public CLI."""
import csv
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

exe, inputs, references, out, tmp = map(Path, sys.argv[1:])
with references.open() as fp:
    refs = {r['case']: r for r in csv.DictReader(fp)}


def get(tag, record):
    path = inputs / f'{tag}.fa'
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != refs[tag]['input_sha256']:
        raise ValueError(tag)
    parts = raw.split(b'>')[1:]
    return b''.join(parts[record].splitlines()[1:])


cases = {}
with (out / 'inputs.csv').open('w') as fp:
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(['case', 'accession_a', 'start_a', 'end_a', 'accession_b', 'start_b',
                     'end_b', 'm', 'n', 'sha256'])
    for size, locus in [(131072, 'insertion'), (262144, 'insertion'), (262144, 'collinear')]:
        left = f'collinear-{size}'
        right = f'{locus}-{size}'
        case = f'{"matched" if locus == "collinear" else "different"}-{size}'
        a, b = get(left, 0), get(right, 1)
        raw = b'>GRCh38\n' + a + b'\n>CHM13\n' + b + b'\n'
        path = tmp / f'{case}.fa'
        path.write_bytes(raw)
        cases[case] = path
        writer.writerow([case, *(refs[left][k] for k in ['accession_a', 'start_a', 'end_a']),
                         *(refs[right][k] for k in ['accession_b', 'start_b', 'end_b']),
                         len(a), len(b), hashlib.sha256(raw).hexdigest()])

rng = random.Random(730519)
scores = {}
with (out / 'loci.csv').open('w') as fp:
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(['block', 'order', 'case', 'algorithm', 'distance', 'seconds', 'peak_rss_kib'])
    for block in range(6):
        configs = [(case, name) for case in cases for name in ('nvpd-cuda-local', 'edlib', 'wfa2')]
        rng.shuffle(configs)
        for order, (case, name) in enumerate(configs, 1):
            record = json.loads(subprocess.check_output(
                [exe, 'pair', '--input', cases[case], '--algorithm', name], text=True))
            score = record['distance']
            if scores.setdefault(case, score) != score:
                raise ValueError(f'{name} disagrees on {case}')
            writer.writerow([block, order, case, name, score, record['seconds'], record['peak_rss_kib']])
            fp.flush()
