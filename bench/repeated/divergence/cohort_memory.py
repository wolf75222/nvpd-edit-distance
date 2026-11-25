#!/usr/bin/env python3
"""Compare one shared CUDA workspace with five simultaneous CUDA processes on one GPU."""
import argparse
from contextlib import ExitStack
import csv
import os
from pathlib import Path
import random
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter_ns

from focus import prepare, sha256


def memory_snapshot(processes):
    rss = pss = 0
    for proc in processes:
        fields = {}
        with Path(f"/proc/{proc.pid}/smaps_rollup").open() as stream:
            for line in stream:
                name, _, value = line.partition(":")
                if name in ("Rss", "Pss"):
                    number, unit = value.split()
                    if unit != "kB":
                        raise ValueError("Unexpected Linux memory unit")
                    fields[name] = int(number)
        if fields["Rss"] < fields["Pss"] or fields["Pss"] <= 0:
            raise ValueError(f"Invalid resident memory for {proc.pid}")
        rss += fields["Rss"]
        pss += fields["Pss"]
    return rss, pss


def dispatch(processes, counts):
    start = perf_counter_ns()
    for proc in processes:
        proc.stdin.write("run\n")
        proc.stdin.flush()
    output = [proc.stdout.readline() for proc in processes]
    seconds = (perf_counter_ns() - start) / 1e9
    peaks = []
    for response, count in zip(output, counts):
        scores, device = map(int, response.split(","))
        if scores != count or device < 0:
            raise RuntimeError(f"Invalid cohort worker reply: {response}")
        peaks.append(device)
    return seconds, sum(peaks)


def batch(worker, mode, cases, directory, cpus):
    groups = [cases] if mode == "shared" else [[case] for case in cases]
    with ExitStack() as stack:
        processes = []
        for cpu, group in zip(cpus, groups):
            args = ["taskset", "-c", str(cpu), worker]
            for item in group:
                args.extend([directory / item["path"], str(item["distance"])])
            proc = stack.enter_context(subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                                       text=True, bufsize=1))
            # Close stdin before Popen.__exit__ waits, including on a protocol error.
            stack.callback(proc.stdin.close)
            processes.append(proc)
        if any(proc.stdout.readline().strip() != "ready" for proc in processes):
            raise RuntimeError("Cohort worker did not load its FASTA")
        counts = [len(group) for group in groups]
        warmup, warm_device = dispatch(processes, counts)
        rss, pss = memory_snapshot(processes)
        duration, device = dispatch(processes, counts)
        for proc in processes:
            proc.stdin.close()
        if any(proc.wait() != 0 for proc in processes):
            raise RuntimeError("Cohort worker failed")
        return {"warmup_seconds": warmup, "seconds": duration, "sum_pss_kib": pss, "sum_rss_kib": rss,
                "warmup_sum_worker_peak_device_bytes": warm_device, "sum_worker_peak_device_bytes": device}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("worker", type=Path)
    parser.add_argument("editdist", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--blocks", type=int, default=5)
    args = parser.parse_args()
    if args.blocks < 2:
        parser.error("Use at least two paired blocks")
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < 6:
        parser.error("Use six Grace cores for the observer and five workers")
    worker = args.worker.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = "1"
    with TemporaryDirectory(prefix="nvpd-cohort-") as temp:
        directory = Path(temp)
        cases = prepare(directory, ["original-five"], args.editdist.resolve())
        os.sched_setaffinity(0, {cpus[-1]})
        ordered_cpus = cpus[:5]
        with (args.output / "cohort.csv").open("x", newline="") as stream:
            writer = None
            rng = random.Random(72988)
            for block in range(1, args.blocks + 1):
                rng.shuffle(cases)
                modes = ("shared", "concurrent") if block % 2 else ("concurrent", "shared")
                for position, mode in enumerate(modes, 1):
                    row = {"job": os.environ.get("SLURM_JOB_ID", "local"), "worker_sha256": sha256(worker),
                           "block": block, "order": position, "mode": mode, "workers": 1 if mode == "shared" else 5,
                           "rep_order": " ".join(r["rep"] for r in cases),
                           "worker_cpus": " ".join(map(str, ordered_cpus[:1] if mode == "shared" else ordered_cpus)),
                           "observer_cpu": cpus[-1], **batch(worker, mode, cases, directory, ordered_cpus)}
                    if writer is None:
                        writer = csv.DictWriter(stream, fieldnames=row, lineterminator="\n")
                        writer.writeheader()
                    writer.writerow(row)
                    stream.flush()
                ordered_cpus = ordered_cpus[1:] + ordered_cpus[:1]
                print(f"Finished cohort block {block}/{args.blocks}", flush=True)


if __name__ == "__main__":
    main()
