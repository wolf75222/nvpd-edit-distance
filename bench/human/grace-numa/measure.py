#!/usr/bin/env python3
"""Measure fixed-work CHM13 alignments with controlled page ownership on romeo-a002."""

import argparse
import csv
import math
import os
from pathlib import Path
import random
import socket
import subprocess

MODES = ("single", "shared0", "shared1", "local", "remote", "twins0", "twins1")
CPUS = (0, 72)
SEED = 840841
TRIALS = ("block", "position", "mode", "workers", "scores", "bases", "warmup_s", "seconds", "cpu_ids")
RESIDENCY = ("block", "mode", "phase", "copy", "chrom", "side", "owner_cpu", "expected_node",
             "observed_node", "bytes", "mapped_bytes", "resident_bytes")


def design():
    """Williams design for seven treatments: balanced positions and directed precedences."""
    rng = random.Random(SEED)
    labels = list(MODES)
    rng.shuffle(labels)
    first = (0, 1, 6, 2, 5, 3, 4)
    orders = [tuple(labels[(i + k) % len(MODES)] for i in first) for k in range(len(MODES))]
    orders += [order[::-1] for order in orders]
    rng.shuffle(orders)
    return orders


def panel_owners(mode):
    return {"single": (0,), "shared0": (0,), "shared1": (72,),
            "local": CPUS, "remote": CPUS[::-1], "twins0": (0, 0), "twins1": (72, 72)}[mode]


def input_lengths(manifest):
    with manifest.open() as stream:
        entries = list(csv.reader(stream, delimiter="\t"))
    if len(entries) != 4 or any(len(row) != 4 for row in entries):
        raise ValueError("Expected four chromosome pairs (FASTA, length, length, distance)")
    return {(i, side): int(row[1 + side]) for i, row in enumerate(entries) for side in (0, 1)}


def topology(output):
    host = socket.gethostname().split(".")[0]
    allowed = sorted(os.sched_getaffinity(0))
    if host != "romeo-a002" or not set(CPUS).issubset(allowed):
        raise RuntimeError(f"Expected CPUs {CPUS} on romeo-a002, got {host} / {allowed}")
    nodes, packages = {}, {}
    with (output / "topology.txt").open("x") as stream:
        print(f"host={host}\naffinity={allowed}\nbase_page_bytes={os.sysconf('SC_PAGE_SIZE')}", file=stream)
        for cpu in CPUS:
            path = Path(f"/sys/devices/system/cpu/cpu{cpu}")
            links = list(path.glob("node[0-9]*"))
            if len(links) != 1:
                raise RuntimeError(f"Cannot resolve NUMA node for CPU {cpu}: {links}")
            nodes[cpu] = int(links[0].name[4:])
            packages[cpu] = int((path / "topology/physical_package_id").read_text())
            core = int((path / "topology/core_id").read_text())
            distances = (links[0] / "distance").read_text().strip()
            print(f"cpu={cpu} core={core} package={packages[cpu]} node={nodes[cpu]} distances={distances}", file=stream)
        print("numa_balancing=" + Path("/proc/sys/kernel/numa_balancing").read_text().strip(), file=stream)
    if len(set(nodes.values())) != 2 or len(set(packages.values())) != 2:
        raise RuntimeError(f"CPUs do not span two Grace chips: {packages}, {nodes}")
    return nodes


def snapshot(procdir, buffers, page_size):
    """Count all resident bytes in the dedicated page-aligned mapping for each sequence."""
    maps = {}
    for line in (procdir / "maps").read_text().splitlines():
        span, flags, *_ = line.split()
        start, end = (int(v, 16) for v in span.split("-"))
        if flags.startswith("r") and flags.endswith("p"):
            maps[start] = end
    counts = {}
    for line in (procdir / "numa_maps").read_text().splitlines():
        address, _, *fields = line.split()
        cells = dict(f.split("=", 1) for f in fields if "=" in f)
        page_kib = int(cells.get("kernelpagesize_kB", 0))
        counts[int(address, 16)] = {int(k[1:]): int(v) * 1024 * page_kib
                                    for k, v in cells.items() if k.startswith("N") and k[1:].isdigit()}
    for copy, chrom, side, cpu, address, length in buffers:
        end = address + ((length + page_size - 1) // page_size) * page_size
        offset = address
        resident = {}
        while offset < end:
            # NUMA policies may split a guarded mapping into adjacent VMAs.
            stop = maps.get(offset)
            if stop is None or stop > end:
                raise RuntimeError(f"Buffer {copy}/{chrom}/{side} lacks a dedicated readable mapping")
            for node, size in counts.get(offset, {}).items():
                resident[node] = resident.get(node, 0) + size
            offset = stop
        mapped = end - address
        if sum(resident.values()) != mapped:
            raise RuntimeError(f"Buffer {copy}/{chrom}/{side} has {sum(resident.values())} resident bytes of {mapped}")
        yield (copy, chrom, side, cpu, length, mapped, resident)


def read_header(child, mode, lengths):
    owners = panel_owners(mode)
    buffers = []
    for line in child.stdout:
        fields = line.split()
        if fields[:1] == ["ready"]:
            break
        if len(fields) != 7 or fields[0] != "buffer":
            raise RuntimeError(f"Unexpected buffer descriptor: {line!r}")
        buffers.append(tuple(int(x, 0) for x in fields[1:]))
    else:
        raise RuntimeError(f"Child exited before warm-up: {child.poll()}")
    expected = {(c, k, s, cpu, length) for c, cpu in enumerate(owners) for (k, s), length in lengths.items()}
    observed = {(c, k, s, cpu, length) for c, k, s, cpu, _, length in buffers}
    if len(buffers) != len(expected) or observed != expected or len({b[4] for b in buffers}) != len(buffers):
        raise RuntimeError("Incomplete, aliased or misplaced panel descriptors")
    workers = 1 if mode == "single" else 2
    bases = 16 * sum(lengths.values())
    if len(fields) != 5 or tuple(map(int, fields[1:4])) != (workers, 64, bases):
        raise RuntimeError(f"Unexpected measured work: {fields}")
    warmup = float(fields[-1])
    if not math.isfinite(warmup) or warmup <= 0:
        raise RuntimeError("Invalid warm-up duration")
    return buffers, (workers, 64, bases, warmup)


def trial(binary, manifest, mode, lengths, nodes, record):
    env = dict(os.environ, OMP_PLACES="{0},{72}", OMP_PROC_BIND="close", OMP_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE", OMP_MAX_ACTIVE_LEVELS="1")
    with subprocess.Popen([binary, manifest, mode, *map(str, CPUS)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          text=True, env=env) as child:
        try:
            buffers, metadata = read_header(child, mode, lengths)
            for phase in ("before", "after"):
                samples = snapshot(Path("/proc") / str(child.pid), buffers, os.sysconf("SC_PAGE_SIZE"))
                for copy, chrom, side, cpu, length, mapped, resident in samples:
                    for node, size in sorted(resident.items()):
                        record((phase, copy, chrom, side, cpu, nodes[cpu], node, length, mapped, size))
                if phase == "before":
                    child.stdin.write("run\n")
                    child.stdin.flush()
                    fields = child.stdout.readline().split()
                    if len(fields) != 2 or fields[0] != "result":
                        raise RuntimeError(f"Missing timing: {fields}")
                    seconds = float(fields[1])
                    if not math.isfinite(seconds) or seconds <= 0:
                        raise RuntimeError("Invalid timed duration")
            child.stdin.write("done\n")
            child.stdin.flush()
            if child.wait() != 0:
                raise RuntimeError("Grace benchmark failed after measurement")
            return (*metadata, seconds, "0" if mode == "single" else "0|72")
        finally:
            if child.poll() is None:
                child.kill()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    nodes = topology(args.output)
    lengths = input_lengths(args.manifest)
    with (args.output / "trials.csv").open("x", newline="") as times, (args.output / "residency.csv").open("x", newline="") as pages:
        clock, residency = csv.writer(times), csv.writer(pages)
        clock.writerow(TRIALS)
        residency.writerow(RESIDENCY)
        for block, order in enumerate(design(), 1):
            for pos, mode in enumerate(order, 1):
                def record(values):
                    residency.writerow((block, mode, *values))
                values = trial(args.binary.resolve(), args.manifest, mode, lengths, nodes, record)
                clock.writerow((block, pos, mode, *values))
                times.flush()
                pages.flush()
                print(f"block={block} position={pos} mode={mode} seconds={values[-2]:.6f}", flush=True)


if __name__ == "__main__":
    main()
