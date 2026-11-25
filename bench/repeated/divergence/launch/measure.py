#!/usr/bin/env python3
"""Compare five score-only implementations on the historical 100-knt pairs."""
import argparse
import csv
import math
import os
from pathlib import Path
import random
import subprocess
import sys
from tempfile import TemporaryDirectory

sys.path.append(str(Path(__file__).resolve().parents[1]))
from cohort_memory import memory_snapshot
from focus import pin_timing_process, prepare, sha256

METHODS = ("nvpd-cuda", "nvpd-cuda-persistent", "nvpd-hybrid-score", "wfa2", "edlib")
GROUPS = ("100k-2pct", "control-5pct", "original-five")


def score(proc, item, method):
    proc.stdin.write("run\n")
    proc.stdin.flush()
    values = proc.stdout.readline().strip().split(",")
    if len(values) != (7 if method == "nvpd-hybrid-score" else 6):
        raise RuntimeError(f"Unexpected worker reply for {method}: {values}")
    m, n, distance, seconds, _, device = values[:6]
    backend = values[6] if method == "nvpd-hybrid-score" else method
    seconds, device = float(seconds), int(device)
    if (int(m) != int(item["length"]) or int(n) != int(item["mutant_length"])
            or int(distance) != item["distance"] or not math.isfinite(seconds) or seconds <= 0
            or device < 0 or (method == "nvpd-hybrid-score"
                              and backend not in ("nvpd-cuda", "nvpd-cuda-persistent", "wfa2"))):
        raise ValueError(f"Invalid result: {item['path']}, {method}, {values}")
    return seconds, device, backend


def measure(binary, method, directory, item):
    args = [binary, directory / item["path"], method, str(item["distance"])]
    if method == "nvpd-hybrid-score":
        args.append("--backend")
    with subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1) as proc:
        if proc.stdout.readline().strip() != "ready":
            raise RuntimeError(f"Worker failed to load {item['path']}")
        warmup, warm_device, warm_backend = score(proc, item, method)
        rss, pss = memory_snapshot([proc])
        seconds, device, backend = score(proc, item, method)
        proc.stdin.close()
        if proc.wait() or warm_backend != backend:
            raise RuntimeError(f"Worker or backend differs between calls for {method}, {item['path']}")
    return dict(backend=backend, warmup_seconds=warmup, seconds=seconds, after_warmup_rss_kib=rss,
                after_warmup_pss_kib=pss, warmup_device_bytes=warm_device, device_bytes=device)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "candidate", "editdist", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--include-extra", action="store_true", help="also run replicates 6 to 10 at 10 percent")
    args = parser.parse_args()
    groups = GROUPS + (("held-out",) if args.include_extra else ())
    binaries = {m: (args.baseline if m == "nvpd-cuda" else args.candidate).resolve() for m in METHODS}
    cli = args.editdist.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.update(OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_PLACES="cores", OMP_PROC_BIND="close")
    pin_timing_process(args.output)
    with (args.output / "binaries.tsv").open("x", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(("job", "program", "sha256"))
        for name, path in (("baseline", args.baseline), ("candidate", args.candidate), ("oracle", cli)):
            w.writerow((os.environ.get("SLURM_JOB_ID", "local"), name, sha256(path)))
    with TemporaryDirectory(prefix="nvpd-launch-") as temp:
        directory = Path(temp)
        items = prepare(directory, groups, cli)
        with (args.output / "inputs.tsv").open("x", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=items[0], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(items)
        rng = random.Random(211785)
        with (args.output / "trials.csv").open("x", newline="") as fh:
            writer = None
            for block in range(len(METHODS)):
                order = list(enumerate(items))
                rng.shuffle(order)
                for index, item in order:
                    # With five cyclic offsets, every method occupies every position once per pair.
                    for position, offset in enumerate((0, 1, 4, 2, 3), 1):
                        method = METHODS[(index + block + offset) % len(METHODS)]
                        row = dict(group=item["group"], rep=item["rep"], block=block+1, order=position,
                                   method=method, distance=item["distance"],
                                   **measure(binaries[method], method, directory, item))
                        if writer is None:
                            writer = csv.DictWriter(fh, fieldnames=row, lineterminator="\n")
                            writer.writeheader()
                        writer.writerow(row)
                    fh.flush()
                print(f"Completed score block {block+1}/5", flush=True)


if __name__ == "__main__":
    main()
