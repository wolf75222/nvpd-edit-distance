#!/usr/bin/env python3
"""Run paired end-to-end CHM13 alignments on Grace, including FASTA decoding and replication."""

import argparse
import csv
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import time
from statistics import mean

# Reuse the previous campaign's checks of the fixed physical CPUs and chromosome lengths.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from measure import CPUS, input_lengths, topology

MODES = ("shared", "local")
PREFIX_BLOCKS = 12
MAX_K = 128
TIMES = ("total_s", "read_s", "copy_s", "setup_s", "align_s")
FIELDS = ("phase", "block", "order", "position", "mode", "rotation", "k", *TIMES, "scores", "process_s")
SEED = 840942


def balanced(blocks, key):
    """Balance the four rotations and AB/BA order, independently of the timings."""
    labels = [(i % 4, (i // 4 + i % 4) % 2) for i in range(blocks)]
    random.Random(SEED + 177 * key + blocks).shuffle(labels)
    return labels


def timing_row(fields, expected_k):
    if len(fields) != 8 or fields[0] != "timing" or int(fields[1]) != expected_k or int(fields[-1]) != expected_k:
        raise ValueError(f"Invalid completed work for K={expected_k}: {fields}")
    t = dict(zip(TIMES, map(float, fields[2:7])))
    if not all(math.isfinite(x) and x >= 0 for x in t.values()) or min(t["total_s"], t["align_s"], t["read_s"]) <= 0:
        raise ValueError(f"Invalid durations: {fields}")
    if (abs(t["total_s"] - t["setup_s"] - t["align_s"]) > 1e-5 or
            t["read_s"] + t["copy_s"] > t["setup_s"] + 1e-5):
        raise ValueError(f"Inconsistent timing boundaries: {fields}")
    return t


def trial(binary, manifest, mode, kind, k, rotation, env):
    cmd = [str(binary), str(manifest), mode, kind, str(k), str(rotation)]
    start = time.perf_counter()
    output = subprocess.check_output(cmd, text=True, env=env)
    process_s = time.perf_counter() - start
    lines = output.splitlines()
    expected = list(range(2, k + 1, 2)) if kind == "prefix" else [k]
    if len(lines) != len(expected):
        raise ValueError(f"Expected {len(expected)} timings, got {len(lines)}")
    results = []
    for n, line in zip(expected, lines):
        t = timing_row(line.split(), n)
        if t["total_s"] > process_s:
            raise ValueError("Native clock exceeds enclosing process duration")
        results.append((n, dict(t, process_s=process_s if kind == "direct" else "")))
    return results


def page_ranges(proc):
    """Include entire VMAs for node counts; pagemap will verify residency on the exact slice."""
    areas = {}
    for row in (proc / "maps").read_text().splitlines():
        span, perms, *_ = row.split()
        if perms.startswith("r"):
            begin, end = (int(x, 16) for x in span.split("-"))
            areas[begin] = end
    counts = {}
    for row in (proc / "numa_maps").read_text().splitlines():
        addr, _, *fields = row.split()
        entries = dict(x.split("=", 1) for x in fields if "=" in x)
        unit = int(entries.get("kernelpagesize_kB", 0)) * 1024
        counts[int(addr, 16)] = {int(key[1:]): int(value) * unit for key, value in entries.items()
                                 if key[0] == "N" and key[1:].isdigit()}
    return areas, counts


def audit_buffers(pid, desc, nodes):
    proc = Path("/proc") / str(pid)
    page = os.sysconf("SC_PAGE_SIZE")
    areas, counts = page_ranges(proc)
    with (proc / "pagemap").open("rb") as pagemap:
        for copy, chrom, side, cpu, address, length in desc:
            start = address // page * page
            stop = ((address + length + page - 1) // page) * page
            pagemap.seek(start // page * 8)
            entries = struct.iter_unpack("=Q", pagemap.read((stop - start) // page * 8))
            resident = sum(bool(x & (1 << 63)) for (x,) in entries) * page
            if resident != stop - start:
                raise ValueError(f"Nonresident chromosome buffer {copy}/{chrom}/{side}")
            foreign = covered = 0
            for base, end in areas.items():
                overlap = min(end, stop) - max(base, start)
                if overlap <= 0:
                    continue
                if not counts.get(base):
                    raise ValueError("No NUMA counts for a resident chromosome VMA")
                # Parsed std::strings may share a larger VMA; attribute all its foreign
                # pages to this buffer for a conservative upper bound, capped at the overlap.
                foreign += min(overlap, sum(n for node, n in counts[base].items() if node != nodes[cpu]))
                covered += overlap
            if covered != resident or foreign >= .01 * resident:
                raise ValueError(f"NUMA ownership unresolved for buffer {copy}/{chrom}/{side}: "
                                 f"{foreign}/{resident} bytes possibly off-node")
            yield copy, chrom, side, cpu, nodes[cpu], length, resident, foreign


def audit(binary, manifest, out, nodes, env):
    lengths = input_lengths(manifest)
    with (out / "residency.csv").open("x", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(("mode", "copy", "chrom", "side", "owner_cpu", "expected_node", "length_bytes",
                         "resident_bytes", "foreign_bytes_upper"))
        for mode in MODES:
            cmd = [str(binary), str(manifest), mode, "direct", "8", "0", "--audit"]
            with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, env=env) as child:
                try:
                    desc, timing = [], False
                    for line in child.stdout:
                        f = line.split()
                        if f == ["audit"]:
                            break
                        if f and f[0] == "timing" and not timing:
                            timing_row(f, 8)
                            timing = True
                        elif len(f) == 7 and f[0] == "buffer":
                            desc.append(tuple(map(int, f[1:])))
                        else:
                            raise ValueError(f"Invalid audit observation: {f}")
                    else:
                        raise ValueError("Native residency handshake missing")
                    copies = 1 if mode == "shared" else 2
                    expected = {(c, i, s, CPUS[c], n) for c in range(copies) for (i, s), n in lengths.items()}
                    if not timing or len(desc) != len(expected) or {(c, i, s, h, n) for c, i, s, h, _, n in desc} != expected:
                        raise ValueError("Incomplete chromosome identities during residency audit")
                    for row in audit_buffers(child.pid, desc, nodes):
                        writer.writerow((mode, *row))
                    child.stdin.write("done\n")
                    child.stdin.flush()
                    if child.wait() != 0:
                        raise ValueError("Native residency audit failed")
                finally:
                    if child.poll() is None:
                        child.kill()


def confirm_plan(rows):
    """Select the direct experiment using independent, even-prefix pilot timings."""
    grouped = {(k, mode): [] for k in range(2, MAX_K + 1, 2) for mode in MODES}
    for row in rows:
        grouped[int(row["k"]), row["mode"]].append(float(row["total_s"]))
    if any(len(x) != PREFIX_BLOCKS for x in grouped.values()):
        raise ValueError("Incomplete calibration")
    means = {k: mean(v) for k, v in grouped.items()}
    even = list(range(2, MAX_K + 1, 2))
    gains = [means[k, "shared"] - means[k, "local"] for k in even]
    cross = next((k for i, k in enumerate(even) if min(gains[i:]) > 0), None)
    center = cross or MAX_K
    keys = set(range(1, 13)) | {32, 64, MAX_K}
    if center <= 32:
        keys.update(range(1, center + 7))
    else:
        keys.update(range(max(1, center - 8), min(MAX_K, center + 6) + 1))
        keys.update(range(16, center, 8))
    keys = sorted(keys)
    # Budget estimate includes the observed startup+score for each fresh process.
    pair_seconds = sum(means[min(MAX_K, k + k % 2), m] for k in keys for m in MODES)
    pilot_seconds = PREFIX_BLOCKS * sum(means[MAX_K, m] for m in MODES)
    blocks = 16 if pilot_seconds + 16 * pair_seconds < 50 * 60 else 12
    return cross, keys, blocks, pilot_seconds + blocks * pair_seconds


def collect(binary, manifest, out, env):
    with (out / "timings.csv").open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        prefix = []
        def pair(kind, length, b, order, assignment):
            rotation, reverse = assignment
            for pos, mode in enumerate(MODES[::-1] if reverse else MODES, 1):
                for k, t in trial(binary, manifest, mode, kind, length, rotation, env):
                    row = dict(phase=kind, block=b, order=order, position=pos, mode=mode,
                               rotation=rotation, k=k, **t, scores=k)
                    writer.writerow(row)
                    if kind == "prefix":
                        prefix.append(row)
                f.flush()
        for b, assignment in enumerate(balanced(PREFIX_BLOCKS, 0), 1):
            pair("prefix", MAX_K, b, 1, assignment)
            print(f"prefix {b}/{PREFIX_BLOCKS}", flush=True)
        cross, keys, blocks, estimate = confirm_plan(prefix)
        with (out / "confirm-design.csv").open("x", newline="") as plan:
            w = csv.writer(plan, lineterminator="\n")
            w.writerow(("pilot_sustained_k", "direct_blocks", "k", "estimated_campaign_s"))
            for k in keys:
                w.writerow((cross, blocks, k, estimate))
        print(f"direct K={keys} blocks={blocks} estimated total {estimate:.0f}s", flush=True)
        assignments = {k: balanced(blocks, k) for k in keys}
        for b in range(1, blocks + 1):
            order = list(keys)
            random.Random(SEED + 1000 + b).shuffle(order)
            for pos, k in enumerate(order, 1):
                pair("direct", k, b, pos, assignments[k][b - 1])
            print(f"direct block {b}/{blocks}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    nodes = topology(args.output)
    env = dict(os.environ, OMP_PLACES="{0},{72}", OMP_PROC_BIND="close", OMP_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE", OMP_MAX_ACTIVE_LEVELS="1")
    binary, manifest = args.binary.resolve(), args.manifest.resolve()
    audit(binary, manifest, args.output, nodes, env)
    collect(binary, manifest, args.output, env)
