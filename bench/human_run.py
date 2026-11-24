#!/usr/bin/env python3
"""Time exact library and CUDA alignments on the pinned chr21 FASTAs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--data", default=Path("data/human"), type=Path)
    p.add_argument("--manifest", default=Path("results/human/manifest.tsv"), type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()
    binary = args.binary.resolve()
    build_hash = digest(binary)
    with args.manifest.open() as manifest:
        cases = list(csv.DictReader(manifest, delimiter="\t"))
    columns = ("job", "case", "algorithm", "trial", "m", "n", "distance", "seconds",
               "wall_seconds", "peak_rss_kib", "device_bytes", "status", "input_sha256", "binary_sha256")
    deadline = time.monotonic() + 32 * 60
    failed = False
    with args.output.open("x", newline="") as target:
        writer = csv.DictWriter(target, columns)
        writer.writeheader()
        for case in cases:
            name = case["case"]
            full = name == "chr21-full-insert64"
            expected = 64 if full else None
            path = args.data / (name + ".fa")
            if digest(path) != case["sha256"]:
                raise ValueError(f"unexpected input digest for {name}")
            backends = ("edlib", "wfa2", "banded-trim") if full else ("edlib", "wfa2", "nvpd-cuda")
            for trial in (1, 2):
                for algorithm in backends:
                    row = dict(job=os.environ.get("SLURM_JOB_ID", "local"), case=name,
                               m=case["m"], n=case["n"], algorithm=algorithm, trial=trial,
                               input_sha256=case["sha256"], binary_sha256=build_hash, status="ok")
                    start = time.monotonic()
                    timeout = min(420, deadline - start)
                    try:
                        if timeout <= 0:
                            raise TimeoutError
                        result = subprocess.run(
                            [str(binary), "pair", "--algorithm", algorithm, "--input", str(path),
                             "--threads", str(args.threads)],
                            text=True, capture_output=True, timeout=timeout, check=True)
                        fields = json.loads(result.stdout)
                        for key in ("distance", "seconds", "peak_rss_kib", "device_bytes"):
                            row[key] = fields[key]
                        if (fields["m"], fields["n"]) != (int(case["m"]), int(case["n"])):
                            raise ValueError(f"length mismatch for {name}")
                        if expected is None:
                            expected = fields["distance"]
                        if fields["distance"] != expected:
                            row["status"] = f"distance differs from {expected}"
                    except (TimeoutError, subprocess.TimeoutExpired):
                        row["status"] = "timeout"
                    except subprocess.CalledProcessError as exc:
                        row["status"] = f"exit {exc.returncode}"
                        print(exc.stderr, flush=True)
                    row["wall_seconds"] = time.monotonic() - start
                    failed |= row["status"] != "ok"
                    writer.writerow(row)
                    target.flush()
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
