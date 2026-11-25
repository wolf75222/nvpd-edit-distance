"""Measure a tighter constructive bound independently of geometry selection."""
import csv
from pathlib import Path
from statistics import mean
import sys

from focus import pin_timing_process

from corridor_sweep import run


def main():
    work, out, large, score_large, control, score_control, plain = sys.argv[1:]
    work, out = Path(work), Path(out)
    binding = out / 'anchor-binding'
    binding.mkdir()
    pin_timing_process(binding)
    names = (plain, plain + '-anchor')
    cases = (('selection', large, score_large), ('control', control, score_control))
    expected = {}
    for name in names:
        with (out / f'plan-{name}.csv').open() as stream:
            for (case, _, _), row in zip(cases, csv.DictReader(stream)):
                if row['policy'] != 'corridor':
                    raise ValueError(f'{name} does not exercise the NvPD corridor on {case}')
                expected[case, name] = int(row['device_bytes'])
    times = {}
    with (out / 'anchor-trials.csv').open('x', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(('case', 'geometry', 'block', 'order', 'score', 'warmup_s', 'time_s', 'device_bytes'))
        for block in range(6):
            for case, path, oracle in cases[::1 if block % 2 == 0 else -1]:
                for pos, name in enumerate(names[::1 if block % 2 == 0 else -1], 1):
                    size = expected[case, name]
                    warm, measured = run(work / name, path, 'nvpd-cuda', oracle, size)
                    times.setdefault((case, name), []).append(measured)
                    writer.writerow((case, name, block + 1, pos, oracle, warm, measured, size))
                    stream.flush()
    ratios = {}
    with (out / 'anchor-selection.csv').open('x', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(('case', 'mean_greedy_s', 'mean_anchored_s', 'ratio'))
        for case, _, _ in cases:
            previous, anchored = (mean(times[case, name]) for name in names)
            ratios[case] = anchored / previous
            writer.writerow((case, previous, anchored, ratios[case]))
    # The historical control may slow down by at most 1% only if the separate
    # 1M selection pair improves by at least 5%.
    use = ((ratios['selection'] <= 1 and ratios['control'] <= 1)
           or (ratios['selection'] <= .95 and ratios['control'] <= 1.01))
    (out/'selected-variant.txt').write_text(names[int(use)]+'\n')
    rows, columns = plain.split('x')
    flags = f'-DCHPS_NVPD_CORRIDOR_ROWS={rows} -DCHPS_NVPD_CORRIDOR_COLUMNS={columns} -DCHPS_NVPD_CORRIDOR_ANCHORS={int(use)}'
    (out/'selected-flags.txt').write_text(flags+'\n')


if __name__ == '__main__':
    main()
