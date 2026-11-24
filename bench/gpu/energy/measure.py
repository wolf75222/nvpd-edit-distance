"""Randomized blocks for the chromosome 21 GPU energy experiment."""
import argparse
import csv
import hashlib
import itertools
import random
import subprocess
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('programs', 'inputs', 'references', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--seconds', type=float, default=2.0)
    opt = p.parse_args()
    cases = [f'{region}-{length}' for region in ('collinear', 'insertion')
             for length in (65536, 524288)]
    with opt.references.open() as fp:
        refs = {r['case']: r for r in csv.DictReader(fp)}
    for case in cases:
        source = opt.inputs / (case + '.fa')
        if hashlib.sha256(source.read_bytes()).hexdigest() != refs[case]['input_sha256']:
            raise ValueError(f'FASTA digest mismatch for {case}')
    fieldnames = ('block,order,case,rows,distance,m,n,cold_s,warm_s,iterations,seconds,'
                  'joules,sensor').split(',')
    rng = random.Random(67038)
    with opt.output.open('w') as out:
        writer = csv.writer(out, lineterminator='\n')
        writer.writerow(fieldnames)
        for block in range(1, 6):
            config = list(itertools.product(cases, (4, 8, 16, 32)))
            rng.shuffle(config)
            for order, (case, rows) in enumerate(config, 1):
                line = subprocess.check_output(
                    [opt.programs / f'trial-{rows}', opt.inputs / f'{case}.fa',
                     refs[case]['distance'], str(opt.seconds)], text=True)
                writer.writerow([block, order, case, rows, *line.strip().split(',')])
                out.flush()


if __name__ == '__main__':
    main()
