"""Select a CUDA tile geometry on an independent sequence pair before confirmation."""
import csv
import math
from pathlib import Path
from random import Random
from statistics import mean
import subprocess
import sys

from focus import pin_timing_process


def run(binary, path, algorithm, oracle, expected):
    lines = subprocess.check_output([binary, path, algorithm, oracle], input='run\nrun\n', text=True).splitlines()
    if len(lines) != 3 or lines[0] != 'ready':
        raise RuntimeError(f'Invalid solver output for {binary}')
    records = [row.split(',') for row in lines[1:]]
    if any(row[2] != oracle or int(row[5]) != expected for row in records):
        raise ValueError(f'{binary}: unexpected score, solver policy or memory')
    warmup, measured = (float(r[3]) for r in records)
    if not all(math.isfinite(v) and v > 0 for v in (warmup, measured)):
        raise ValueError(f"{binary}: nonpositive or nonfinite measurement")
    return warmup, measured


def main():
    work, out, large, score_large, control, score_control, *geometries = sys.argv[1:]
    work, out = Path(work), Path(out)
    pin_timing_process(out)
    reference = '32x128'
    expected_bytes = {}
    eligible = set(geometries)
    for geometry in geometries:
        with (out / f'plan-{geometry}.csv').open() as stream:
            for case, row in zip(('selection', 'control'), csv.DictReader(stream)):
                if row['policy'] != 'corridor':
                    eligible.discard(geometry)
                else:
                    expected_bytes[case, geometry] = int(row['device_bytes'])
    if reference not in eligible:
        raise ValueError('Reference geometry does not exercise the corridor')
    order = [g for g in geometries if g in eligible] + ['wfa2']
    rng = Random(17312)
    rng.shuffle(order)
    cases = [('selection', large, score_large), ('control', control, score_control)]
    timings = {}
    with (out/'trials.csv').open('x', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(('case', 'geometry', 'block', 'order', 'score', 'warmup_s', 'time_s', 'device_bytes'))
        for block in range(len(order)):
            rng.shuffle(cases)
            for case, path, oracle in cases:
                for position in range(len(order)):
                    name = order[(position+block) % len(order)]
                    binary = work / (reference if name == 'wfa2' else name)
                    algorithm = 'wfa2' if name == 'wfa2' else 'nvpd-cuda'
                    expected = 0 if name == 'wfa2' else expected_bytes[case, name]
                    warmup, measured = run(binary, path, algorithm, oracle, expected)
                    timings.setdefault((case, name), []).append(measured)
                    writer.writerow((case, name, block+1, position+1, oracle, warmup, measured, expected))
                    stream.flush()
    # Among candidates no slower on the historical control, minimize mean
    # elapsed time on the separate 1M selection pair.
    control_limit = mean(timings['control', reference])
    names = [g for g in geometries if g in eligible and mean(timings['control', g]) <= control_limit]
    winner = min(names, key=lambda g: mean(timings['selection', g]))
    (out/'selected.txt').write_text(winner+'\n')


if __name__ == '__main__':
    main()
