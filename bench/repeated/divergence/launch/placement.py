#!/usr/bin/env python3
"""Measure four fixed chromosome/worker assignments and sample Linux NUMA residency."""
import csv
import os
from pathlib import Path
import random
import subprocess
import sys


def residency(pid, buffers):
    """Count resident pages once per VMA, including VMA slack when it contains a sequence."""
    proc = Path("/proc") / str(pid)
    ranges = {}
    for row in (proc / "maps").read_text().splitlines():
        start, end = (int(value, 16) for value in row.split(maxsplit=1)[0].split("-"))
        ranges[start] = end
    grouped = {}
    found = set()
    for row in (proc / "numa_maps").read_text().splitlines():
        address, _, *fields = row.split()
        start = int(address, 16)
        end = ranges.get(start)
        values = dict(part.split("=", 1) for part in fields if "=" in part)
        if not end or "file" in values or "anon" not in values:
            continue
        page_kib = int(values["kernelpagesize_kB"])
        scopes = []
        if int(values["anon"]) * page_kib >= 16384:
            scopes.append("anonymous_ge_16_MiB")
        overlaps = [i for i, (lo, hi) in enumerate(buffers) if lo < end and start < hi]
        if overlaps:
            scopes.append("sequence_vmas")
            found.update(overlaps)
        for scope in scopes:
            for key, value in values.items():
                if key.startswith("N") and key[1:].isdigit():
                    group = scope, int(key[1:]), page_kib
                    grouped[group] = grouped.get(group, 0) + int(value)
    if len(found) != len(buffers):
        raise ValueError("Could not identify a VMA for every resident sequence")
    return grouped


def run(binary, manifest, out, cpus, env):
    first = ["owner", "master"] * 3
    random.Random(840841).shuffle(first)
    with (out / "trials.csv").open("w", newline="") as trials, (out / "residency.csv").open("w", newline="") as pages:
        times, numa = csv.writer(trials), csv.writer(pages)
        times.writerow(("block", "position", "mode", "pairs", "bases", "read_seconds", "warmup_seconds",
                        "seconds", "cpu_ids"))
        numa.writerow(("block", "mode", "scope", "numa_node", "page_kib", "resident_pages"))
        for block, initial in enumerate(first, 1):
            for position, mode in enumerate((initial, "owner" if initial == "master" else "master"), 1):
                with subprocess.Popen([binary, manifest, mode], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      text=True, env=env) as child:
                    try:
                        row = child.stdout.readline().split()
                        if len(row) != 23 or row[0] != "ready" or [int(v) for v in row[3:7]] != cpus:
                            raise ValueError(f"Invalid worker readiness or affinity: {row[:7]}")
                        read, warmup = map(float, row[1:3])
                        buffers = [(int(row[i]), int(row[i]) + int(row[i+1])) for i in range(7, len(row), 2)]
                        try:
                            snapshot = residency(child.pid, buffers)
                        except (PermissionError, FileNotFoundError) as error:
                            (out / "residency-unavailable.txt").write_text(str(error) + "\n")
                        else:
                            for key, value in sorted(snapshot.items()):
                                numa.writerow((block, mode, *key, value))
                        pages.flush()
                        child.stdin.write("run\n")
                        child.stdin.flush()
                        seconds = float(child.stdout.readline())
                        child.stdin.close()
                        if child.wait() or min(read, warmup, seconds) <= 0:
                            raise ValueError("Incomplete placement trial")
                        bases = 16 * sum(hi - lo for lo, hi in buffers)
                        times.writerow((block, position, mode, 64, bases, read, warmup, seconds,
                                        "|".join(map(str, cpus))))
                        trials.flush()
                    finally:
                        if child.poll() is None:
                            child.kill()


def main():
    binary, previous, out = map(Path, sys.argv[1:])
    cores = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        package, core = (int((topology / name).read_text()) for name in ("physical_package_id", "core_id"))
        cores.setdefault((package, core), cpu)
    # First choose one core per available package, then fill unused places.
    first = {}
    for (package, _), cpu in cores.items():
        first.setdefault(package, cpu)
    cpus = list(dict.fromkeys([*first.values(), *cores.values()]))[:4]
    if len(cpus) != 4:
        raise ValueError(f"Expected four physical cores, available: {cores}")
    env = dict(os.environ, OMP_PLACES=",".join(f"{{{cpu}}}" for cpu in cpus),
               OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close", OMP_MAX_ACTIVE_LEVELS="1")
    with (out / "environment.txt").open("a") as stream:
        print("affinity=" + ",".join(map(str, sorted(os.sched_getaffinity(0)))), file=stream)
        for (package, core), cpu in cores.items():
            if cpu in cpus:
                directory = Path(f"/sys/devices/system/cpu/cpu{cpu}")
                nodes = sorted(path.name for path in directory.glob("node[0-9]*"))
                print(f"cpu={cpu} package={package} core={core} numa={','.join(nodes)}", file=stream)
    run(binary, previous / "inputs.tsv", out, cpus, env)


if __name__ == "__main__":
    main()
