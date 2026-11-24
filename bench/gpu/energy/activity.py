"""Sum CUDA activity on one chromosome window, keeping transfer and kernel time distinct."""
import csv
import sqlite3
import sys
from pathlib import Path

with Path(sys.argv[2]).open('w') as fp:
    writer = csv.writer(fp, lineterminator='\n')
    writer.writerow(('case', 'rows', 'launches', 'kernel_s', 'copies', 'copy_bytes', 'copy_s'))
    for rows in (4, 8, 16, 32):
        with sqlite3.connect(Path(sys.argv[1]) / f'activity-{rows}.sqlite') as con:
            kernels = con.execute('SELECT count(*), sum(end-start)/1e9 '
                                  'FROM CUPTI_ACTIVITY_KIND_KERNEL').fetchone()
            copies = con.execute('SELECT count(*), sum(bytes), sum(end-start)/1e9 '
                                 'FROM CUPTI_ACTIVITY_KIND_MEMCPY').fetchone()
        writer.writerow(('collinear-65536', rows, *kernels, *copies))
