"""Exercise short final tiles and fallback paths, then pair two kernels within each block."""
import csv
import hashlib
import random
import subprocess
import sys
from pathlib import Path


def distance(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        row = [i]
        for j, y in enumerate(b, 1):
            row.append(min(previous[j] + 1, row[-1] + 1, previous[j - 1] + (x != y)))
        previous = row
    return previous[-1]


def main():
    tmp, inputs, references, out = map(Path, sys.argv[1:])
    rng = random.Random(843902)
    programs = ('direct', 'packed')
    tests = []
    for length in (15, 16, 17, 31, 32, 33):
        tests.append(('ACGT' * (length // 4) + 'AGC'[:length % 4],
                      ''.join(rng.choices('ACGT', k=880)) + 'T' * 265 + 'ACG'))
    for gap in (254, 255, 256, 4700, 9000):
        tests.append(('ACGT' * 10 + 'A', 'ACGT' + 'A' * gap + 'TGCA'))
    tests.extend([('ACGTN' * 7, 'ACGTN' + 'G' * 340 + 'CAT'),
                  ('A' * 34, 'A' * 263), ('CGTA', 'ACGT' * 33)])
    with (out / 'validation.csv').open('w') as fp:
        writer = csv.writer(fp, lineterminator='\n')
        writer.writerow(['pair', 'm', 'n', 'distance', *programs])
        for i, (a, b) in enumerate(tests):
            expected = distance(a, b)
            (tmp / 'test.fa').write_text(f'>a\n{a}\n>b\n{b}\n')
            for name in programs:
                subprocess.run([tmp / name, tmp / 'test.fa', str(expected), '0'], check=True)
            writer.writerow([i, len(a), len(b), expected, expected, expected])
    with references.open() as fp:
        refs = {r['case']: r for r in csv.DictReader(fp)}
    cases = [f'{region}-{length}' for region in ('collinear', 'insertion')
             for length in (65536, 524288)]
    for case in cases:
        if hashlib.sha256((inputs / f'{case}.fa').read_bytes()).hexdigest() != refs[case]['input_sha256']:
            raise ValueError(f'incorrect input for {case}')
    with (out / 'packing.csv').open('w') as fp:
        writer = csv.writer(fp, lineterminator='\n')
        writer.writerow(('block,order,case,kernel,distance,m,n,cold_s,warm_s,iterations,'
                         'seconds,joules,sensor').split(','))
        for block in range(1, 6):
            configs = [(case, p) for case in cases for p in programs]
            rng.shuffle(configs)
            for order, (case, p) in enumerate(configs, 1):
                record = subprocess.check_output([tmp / p, inputs / f'{case}.fa',
                                                 refs[case]['distance'], '2'], text=True)
                writer.writerow([block, order, case, p, *record.strip().split(',')])
                fp.flush()


if __name__ == '__main__':
    main()
