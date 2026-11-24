#!/usr/bin/env python3
"""Run one alignment per process so peak RSS belongs to one pair and backend."""
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import socket
import subprocess
import sys
import time
from pathlib import Path

CPU = ("seq", "diagonal", "nvpd-omp", "edlib", "wfa2")
PARALLEL = {"diagonal", "nvpd-omp"}
FIELDS = ("job", "host", "arch", "target", "binary_sha256", "path", "input_sha256",
          "length", "rate", "rep", "seed", "m", "n", "algorithm", "threads", "trial",
          "distance", "seconds", "wall_seconds", "peak_rss_kib", "device_bytes", "status", "detail")


def digest(path):
    with path.open("rb") as stream:
        h = hashlib.sha256()
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--manifest", default=Path("data/synthetic/manifest.tsv"), type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--target", choices=("cpu", "gpu"), required=True)
    p.add_argument("--profile", choices=("quick", "full"), default="quick")
    p.add_argument("--algorithms", nargs="+", choices=(*CPU, "nvpd-cuda"))
    p.add_argument("--threads", type=int, nargs="+", help="thread counts for OpenMP backends")
    p.add_argument("--lengths", type=int, nargs="+")
    p.add_argument("--rates", type=float, nargs="+")
    p.add_argument("--replicates", type=int, nargs="+", help="one-based pair numbers per length/rate")
    p.add_argument("--trials", type=int, default=3)
    p.add_argument("--timeout", type=float, help="seconds per process")
    a = p.parse_args()
    cores = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    threads = a.threads or sorted({1, min(2, cores), min(4, cores), min(8, cores), cores})
    if min(threads) < 1 or max(threads) > cores or a.trials < 1:
        p.error("threads must fit the allocation and trials must be positive")
    if a.timeout is not None and a.timeout <= 0:
        p.error("timeout must be positive")
    algorithms = a.algorithms or (CPU if a.target == "cpu" else ("nvpd-cuda",))
    if "nvpd-cuda" in algorithms and a.target != "gpu":
        p.error("nvpd-cuda needs a GPU target")

    with a.manifest.open(newline="") as stream:
        pairs = list(csv.DictReader(stream, delimiter="\t"))
    lengths = a.lengths or ([5000, 20000] if a.profile == "quick" else None)
    replicates = a.replicates or ([1] if a.profile == "quick" else None)
    pairs = [x for x in pairs if (not lengths or int(x["length"]) in lengths)
             and (not a.rates or float(x["rate"]) in a.rates)
             and (not replicates or int(x["rep"]) in replicates)]
    if algorithms == ["seq"]:
        pairs = [x for x in pairs if int(x["length"]) < 1000000 or int(x["rep"]) == 1]
    if not pairs:
        p.error("no eligible measurements match these filters (1M sequential uses replicate 1)")

    binary = a.binary.resolve(strict=True)
    common = dict(job=os.environ.get("SLURM_JOB_ID", "local"), host=socket.gethostname(),
                  arch=platform.machine(), target=a.target, binary_sha256=digest(binary))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with a.output.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for pair in pairs:
            path = a.manifest.parent / pair["path"]
            meta = {k: pair[k] for k in ("path", "length", "rate", "rep", "seed")}
            meta.update(input_sha256=digest(path), m=int(pair["length"]), n=int(pair["mutant_length"]))
            for algorithm in algorithms:
                large_seq = algorithm == "seq" and int(pair["length"]) >= 1000000
                if large_seq and int(pair["rep"]) != 1:
                    continue
                for workers in (threads if algorithm in PARALLEL else [1]):
                    cmd = [str(binary), "pair", "--input", str(path), "--algorithm", algorithm,
                           "--threads", str(workers)]
                    env = dict(os.environ, OMP_NUM_THREADS=str(workers), OMP_DYNAMIC="FALSE")
                    for trial in range(1, 1 + (1 if large_seq else a.trials)):
                        row = dict(common, **meta, algorithm=algorithm, threads=workers, trial=trial,
                                   status="ok")
                        start = time.monotonic()
                        try:
                            result = subprocess.run(cmd, capture_output=True, text=True, env=env,
                                                    check=True, timeout=a.timeout)
                            stats = json.loads(result.stdout)
                            if (stats["m"], stats["n"]) != (meta["m"], meta["n"]):
                                raise ValueError("FASTA lengths differ from manifest")
                            row.update({k: stats[k] for k in
                                        ("distance", "seconds", "peak_rss_kib", "device_bytes")})
                            if not math.isfinite(row["seconds"]) or row["seconds"] <= 0:
                                raise ValueError("invalid alignment time")
                        except (subprocess.SubprocessError, ValueError, KeyError, OSError) as exc:
                            row["status"] = "timeout" if isinstance(exc, subprocess.TimeoutExpired) else "error"
                            row["detail"] = (getattr(exc, "stderr", None) or str(exc))[-1000:]
                            failed = True
                            print(f"{path} {algorithm}/{workers}: {row['detail']}", file=sys.stderr)
                        row["wall_seconds"] = time.monotonic() - start
                        writer.writerow(row)
                        stream.flush()
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
