"""Recreate plans for measured inputs and verify both the FASTA hash and the recorded GPU allocation."""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile

from focus import sha256
from synthetic import make_pair, write_pair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path, help='compiled corridor_plan executable')
    parser.add_argument('timing', type=Path, help='directory containing inputs.tsv and trials.csv')
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    with (args.timing / 'inputs.tsv').open() as stream:
        inputs = list(csv.DictReader(stream, delimiter='\t'))
    device_bytes, blocks = {}, {}
    with (args.timing / 'trials.csv').open() as stream:
        for row in csv.DictReader(stream):
            if row['method'] == 'candidate':
                key = row['group'], row['rep']
                device_bytes.setdefault(key, set()).add(int(row['device_bytes']))
                observed = blocks.setdefault(key, set())
                if int(row['block']) in observed:
                    raise ValueError(f'Repeated block for {key}')
                observed.add(int(row['block']))
    with tempfile.TemporaryDirectory() as work:
        paths = []
        for item in inputs:
            target = Path(work) / Path(item['path']).name
            edits = sum(int(item[k]) for k in ('insertions', 'deletions', 'substitutions'))
            a, b, _ = make_pair(int(item['length']), edits, int(item['seed']))
            write_pair(target, target.stem, a, b)
            if sha256(target) != item['sha256']:
                raise ValueError(f'FASTA differs from measured input: {target.name}')
            paths.append(target)
        lines = subprocess.check_output([args.plan.resolve(), *paths], text=True).splitlines()
    plans = list(csv.DictReader(lines))
    if len(plans) != len(inputs):
        raise ValueError('Missing corridor plans')
    rows = []
    for item, plan in zip(inputs, plans):
        if plan['path'] != Path(item['path']).name or plan['policy'] != 'corridor':
            raise ValueError(f'Unexpected plan for {item["path"]}')
        key = item['group'], item['rep']
        if blocks.get(key) != set(range(1, 9)) or int(plan['upper']) < int(item['distance']):
            raise ValueError(f'Incomplete measurements or invalid upper bound for {item["path"]}')
        if device_bytes[key] != {int(plan['device_bytes'])}:
            raise ValueError(f'Analytical and measured allocation differ for {item["path"]}')
        rows.append({k: item[k] for k in ('group', 'rep', 'sha256', 'distance')} | plan)
    with args.output.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
