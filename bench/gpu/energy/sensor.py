"""Measure temporal integration sensitivity and inspect the selected GPU paths."""
import csv
import random
import subprocess
import sys
from pathlib import Path

tmp, inputs, out = map(Path, sys.argv[1:])
with (out / 'geometry.csv').open('w') as fp:
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow('case,sms,shared_limit,columns,pitch,reach,index_entries,sentinels,kernel'.split(','))
    for size in (65536, 524288):
        for region in ('collinear', 'insertion'):
            name = f'{region}-{size}'
            geom = subprocess.check_output([tmp / 'geometry', inputs / f'{name}.fa'], text=True)
            writer.writerow([name, *geom.strip().split(',')])
    for gap in (255, 4700, 9000):
        path = tmp / 'fallback.fa'
        path.write_text(f'>a\n{"ACGT" * 10}A\n>b\nACGT{"A" * gap}TGCA\n')
        geom = subprocess.check_output([tmp / 'geometry', path], text=True)
        writer.writerow([f'gap-{gap}', *geom.strip().split(',')])
rng = random.Random(810297)
with (out / 'window.csv').open('w') as fp:
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(('block,order,target_s,distance,m,n,cold_s,warm_s,iterations,'
                     'seconds,joules,sensor').split(','))
    for block in range(1, 6):
        durations = [.5, 2., 8.]
        rng.shuffle(durations)
        for order, sec in enumerate(durations, 1):
            line = subprocess.check_output([tmp / 'trial', inputs / 'collinear-65536.fa',
                                           '51', str(sec)], text=True)
            writer.writerow([block, order, sec, *line.strip().split(',')])
            fp.flush()
