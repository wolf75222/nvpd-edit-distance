#!/usr/bin/env python3
"""Measure complete chromosomes with one fixed physical CPU per lane."""
from contextlib import ExitStack
import argparse
import csv
import hashlib
import os
from pathlib import Path
import random
import subprocess

READ_CASES = ("chr1", "chr4", "chr7", "chr12", "chr16", "chr20")
EVENTS = ("instructions", "cycles", "L1-dcache-loads", "L1-dcache-load-misses")


def run_trial(cmd, case):
    lines = subprocess.check_output(cmd, text=True, timeout=480).splitlines()
    if len(lines) != 2:
        raise ValueError("expected one warm and one timed call")
    parsed = [line.split(",") for line in lines]
    if any((int(row[0]), int(row[1])) != (int(case["m"]), int(case["n"])) for row in parsed):
        raise ValueError("length mismatch")
    if parsed[0][2] != parsed[1][2]:
        raise ValueError("score mismatch")
    return dict(distance=int(parsed[1][2]), warm_seconds=float(parsed[0][3]),
                seconds=float(parsed[1][3]), peak_rss_kib=int(parsed[1][4]))


def writer(stream, fields):
    sink = csv.DictWriter(stream, fields.split(), lineterminator="\n")
    sink.writeheader()
    return sink


def normalized_digest(binary, path):
    with subprocess.Popen([binary, path, "check"], stdout=subprocess.PIPE) as command:
        command.stdout.readline()
        digest = hashlib.sha256()
        for part in iter(lambda: command.stdout.read(1 << 20), b""):
            digest.update(part)
        if command.wait():
            raise ValueError(f"failed to decode {path}")
    return digest.hexdigest()


def canonical_file_digest(path):
    digest = hashlib.sha256()
    first = True
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                if not first:
                    digest.update(b"\n")
                digest.update(line)
                first = False
            else:
                digest.update(line.rstrip(b"\r\n"))
    digest.update(b"\n")
    return digest.hexdigest()


def profile(binary, path):
    rows = subprocess.check_output([binary, path, "profile"], text=True, timeout=240).splitlines()
    counters = {}
    for record in rows[2:]:
        name, count, pct = record.split(",")
        counters[name] = (int(count), float(pct))
    if set(EVENTS) != counters.keys() or any(pct < 99 for _, pct in counters.values()):
        raise ValueError(f"missing or multiplexed events: {counters}")
    return counters


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data", type=Path)
    p.add_argument("build", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    lane = int(os.environ["SLURM_PROCID"])
    lanes = int(os.environ["SLURM_NTASKS"])
    cores = sorted(os.sched_getaffinity(0))
    if len(cores) != 1:
        raise RuntimeError(f"expected a single CPU, got {cores}")
    with (args.data / "manifest.tsv").open() as source:
        cases = list(csv.DictReader(source, delimiter="\t"))
    windows = cases[lane::lanes]
    rng = random.Random(80271 + lane)
    binaries = {name: args.build / ("panel_" + name) for name in ("initial", "scalar", "simd")}
    args.out.mkdir(parents=True, exist_ok=True)
    common = dict(job=os.environ["SLURM_JOB_ID"], lane=lane, cpu=cores[0])
    with ExitStack() as stack:
        times, reads, counters, checks = [stack.enter_context((args.out / f"{kind}-{lane}.csv").open("x"))
                                          for kind in ("trials", "reader", "counters", "checks")]
        keys = "job lane cpu case block method"
        fields = keys + " distance warm_seconds seconds peak_rss_kib"
        time_out, read_out = writer(times, fields), writer(reads, fields)
        count_out = writer(counters, keys + " event value running_pct")
        check_out = writer(checks, "case method normalized_sha256")
        for case in windows:
            name = case["case"]
            path = args.data / (name + ".fa")
            control = None
            for block in range(1, 6):
                methods = ["edlib", "wfa2"]
                rng.shuffle(methods)
                for method in methods:
                    row = run_trial([binaries["simd"], path, method], case)
                    if control is not None and row["distance"] != control:
                        raise ValueError(f"discordant score on {name}")
                    control = row["distance"]
                    time_out.writerow(dict(common, case=name, block=block, method=method, **row))
                    times.flush()
            if name not in READ_CASES:
                continue
            # Check every returned byte; file-PMU runs subsequently use the page cache.
            expected = canonical_file_digest(path)
            for method, binary in binaries.items():
                value = normalized_digest(binary, path)
                if value != expected:
                    raise ValueError(f"decoded chromosome mismatch on {name}, {method}")
                check_out.writerow(dict(case=name, method=method, normalized_sha256=value))
            for block in range(1, 6):
                methods = list(binaries)
                rng.shuffle(methods)
                for method in methods:
                    binary = binaries[method]
                    key = dict(common, case=name, block=block, method=method)
                    read_out.writerow(dict(key, **run_trial([binary, path, "reader"], case)))
                    reads.flush()
                    for event, (value, pct) in profile(binary, path).items():
                        count_out.writerow(dict(key, event=event, value=value, running_pct=pct))
                    counters.flush()


if __name__ == "__main__":
    main()
