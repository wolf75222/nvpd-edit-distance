#!/usr/bin/env python3
"""Profile a cost-only NvPD call and measure strong, weak and paired A/B trials."""
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PHASE, BUILD, OUT = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
JOB = os.environ.get("SLURM_JOB_ID", "local")
THREADS = (1, 2, 4, 8, 16)
PROFILE = (("square", 20000, 20000), ("wide", 256, 1000000), ("over_l3", 64, 4194304))
FIELDS = ["job", "study", "case", "variant", "threads", "trial", "m", "n", "cells",
          "input_fnv64", "reference_distance", "distance", "seconds", "peak_rss_kib",
          "dp_bytes", "mit_bytes", "workspace_peak_bytes", "cpu_ids", "max_gap", "comparison", "source_sha256"]
EVENT_FIELDS = ["job", "case", "variant", "threads", "trial", "m", "n", "input_fnv64",
                "reference_distance", "distance", "event", "value", "unit", "running_pct", "scope", "source_sha256"]
ORACLES = {}


def source_sha(variant):
    return (OUT / ("source-" + variant + ".sha256")).read_text().strip()



def writer(path, fields):
    stream = path.open("x", newline="")
    output = csv.DictWriter(stream, fieldnames=fields)
    output.writeheader()
    return stream, output


def command(variant, m, n, p, warmups=1, algorithm="nvpd"):
    return [str(BUILD / variant), algorithm, str(m), str(n), str(p), str(warmups)]


def run(cmd, env=None):
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=120)
    if proc.returncode:
        raise RuntimeError(f"{Path(cmd[0]).name} exited {proc.returncode}: {proc.stderr[-700:]}")
    fields = proc.stdout.strip().split(",")
    if len(fields) != 6:
        raise RuntimeError(f"missing kernel result: {proc.stderr[-700:]}")
    return fields


def expected(m, n):
    if (m, n) not in ORACLES:
        d, _, _, h, _, _ = run(command("baseline", m, n, 1, 0, "seq"))
        ORACLES[m, n] = int(d), h
    return ORACLES[m, n]


def checked(result, m, n):
    d, t, rss, h, cpus, gap = result
    reference, digest = expected(m, n)
    if (int(d), h) != (reference, digest):
        raise RuntimeError(f"distance/input mismatch at {m}x{n}")
    return dict(m=m, n=n, input_fnv64=h, distance=int(d), reference_distance=reference,
                seconds=float(t), peak_rss_kib=int(rss), cpu_ids=cpus, max_gap=int(gap))


def timings(study, configs):
    stream, out = writer(OUT / (study + ".csv"), FIELDS)
    rng = random.Random(20261124)
    with stream:
        for trial in range(1, 4):
            order = list(configs)
            rng.shuffle(order)
            for case, m, n, p, variants in order:
                # Rotate A/B ordering by block, using a fresh equally warmed process for each.
                for v in (variants if trial % 2 else variants[::-1]):
                    row = checked(run(command(v, m, n, p)), m, n)
                    dp = 8 * (n + 1)
                    integer_index = 16 * (n + 1)
                    offsets = 4 * (n + 1)
                    mit = offsets if v == "v2" and row["max_gap"] <= 255 else integer_index
                    peak = max(dp + mit, integer_index + offsets) if v == "v2" else dp + mit
                    out.writerow(dict(row, job=JOB, study=study if study == "ab" else case,
                                      case=case, variant=v, threads=p, trial=trial, cells=m*n,
                                      dp_bytes=dp, mit_bytes=mit, workspace_peak_bytes=peak,
                                      comparison="-".join(variants), source_sha256=source_sha(v)))
                    stream.flush()


def profiles(variant):
    perf = os.environ.get("PERF_BIN") or shutil.which("perf")
    if not perf:
        raise RuntimeError("Set PERF_BIN to a perf executable supporting stat --control.")
    stream, out = writer(OUT / ("counters-" + variant + ".csv"), EVENT_FIELDS)
    groups = (("cycles", "instructions", "branches", "branch-misses"),
              ("cache-references", "cache-misses", "L1-dcache-loads", "L1-dcache-load-misses"),
              ("dTLB-loads", "dTLB-load-misses"))
    with stream, tempfile.TemporaryDirectory(dir=BUILD) as work:
        wd = Path(work)
        os.mkfifo(wd / "ctl")
        os.mkfifo(wd / "ack")
        ctl = f"fifo:{wd / 'ctl'},{wd / 'ack'}"
        env = dict(os.environ, CHPS_PERF_CTL=str(wd / "ctl"), CHPS_PERF_ACK=str(wd / "ack"), LC_ALL="C")
        perf_version = subprocess.check_output([perf, "--version"], text=True).strip()
        (OUT / "perf-version.txt").write_text(perf_version + "\n")
        for case, m, n in PROFILE:
            for p in (1, 8):
                for trial in (1, 2):
                    for events in groups:
                        pmu = ",".join(x + ":u" for x in events)
                        cmd = [perf, "stat", "-x", ";", "--no-big-num", "--delay=-1",
                               "--control", ctl, "-o", str(wd / "perf.csv"), "-e", pmu,
                               "--", *command(variant, m, n, p)]
                        result = checked(run(cmd, env), m, n)
                        for line in (wd / "perf.csv").read_text().splitlines():
                            cols = line.split(";")
                            if len(cols) < 5 or not cols[2].strip():
                                continue
                            value, unit, event, _, running = [x.strip() for x in cols[:5]]
                            out.writerow(dict(job=JOB, case=case, variant=variant, threads=p, trial=trial,
                                              m=m, n=n, input_fnv64=result["input_fnv64"],
                                              reference_distance=result["reference_distance"], distance=result["distance"],
                                              event=event, value=value, unit=unit or "count", running_pct=running,
                                              scope="warmed_nvpd_call_userspace_all_threads",
                                              source_sha256=source_sha(variant)))
                        stream.flush()


if PHASE == "baseline":
    profiles("baseline")
elif PHASE in ("v1", "v2"):
    variants = ("baseline", "v1") if PHASE == "v1" else ("v1", "v2")
    timings("ab", [(case, m, n, p, variants) for case, m, n in PROFILE for p in (1, 8, 16)])
    profiles(PHASE)
    configs = [("strong", 20000, 20000, p, (PHASE,)) for p in THREADS]
    # Every worker gets exactly 4096*4096 cells and every team crosses 4096 row barriers.
    configs += [("weak", 4096, 4096*p, p, (PHASE,)) for p in THREADS]
    timings("scaling", configs)
else:
    raise SystemExit("phase must be baseline, v1, or v2")
