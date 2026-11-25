#!/usr/bin/env python3
"""Time two NvPD builds, Edlib and WFA2 in paired blocks on synthetic inputs."""
import argparse
import csv
import ctypes
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/datasets"))
from synthetic import generate

METHODS = ("original", "candidate", "edlib", "wfa2")
DEFAULT_GROUPS = ("original-five", "held-out", "control-5pct")
GROUPS = {
    "original-five": (100000, "0.1", range(1, 6)),
    "held-out": (100000, "0.1", range(6, 11)),
    "control-5pct": (100000, "0.05", range(1, 6)),
    "100k-2pct": (100000, "0.02", range(1, 6)),
    "1m-2pct": (1000000, "0.02", range(1, 6)),
}
TRIAL_FIELDS = ("group", "rep", "block", "order", "method", "distance", "warmup_seconds",
                "seconds", "rss_kib", "device_bytes")
MASTER_SEED = 2026


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gpu_pci_address():
    # The driver resolves CUDA_VISIBLE_DEVICES for ordinal zero without a CUDA context.
    try:
        driver = ctypes.CDLL("libcuda.so.1")
        device = ctypes.c_int()
        address = ctypes.create_string_buffer(32)
        if (driver.cuInit(0) or driver.cuDeviceGet(ctypes.byref(device), 0)
                or driver.cuDeviceGetPCIBusId(address, len(address), device)):
            return ""
        return address.value.decode("ascii").lower()
    except (OSError, AttributeError):
        return ""


def pin_timing_process(directory):
    allowed = sorted(os.sched_getaffinity(0))
    pci = gpu_pci_address()
    node = -1
    method = "first-allowed"
    cpu = allowed[0]
    if pci:
        try:
            node = int((Path("/sys/bus/pci/devices") / pci / "numa_node").read_text())
            if node >= 0:
                sysnodes = Path("/sys/devices/system/node")
                distances = [int(s) for s in (sysnodes / f"node{node}/distance").read_text().split()]
                node_ids = sorted(int(p.name[4:]) for p in sysnodes.glob("node[0-9]*"))
                if len(distances) == len(node_ids):
                    costs = dict(zip(node_ids, distances))
                    candidates = []
                    for c in allowed:
                        links = list(Path(f"/sys/devices/system/cpu/cpu{c}").glob("node[0-9]*"))
                        if len(links) == 1:
                            cpu_node = int(links[0].name[4:])
                            candidates.append((costs[cpu_node], c, cpu_node))
                    if candidates:
                        _, cpu, cpu_node = min(candidates)
                        method = "gpu-numa" if node == cpu_node else "nearest-allowed-numa"
        except (OSError, ValueError, KeyError):
            pass
    os.sched_setaffinity(0, {cpu})
    if os.sched_getaffinity(0) != {cpu}:
        raise RuntimeError("Failed to bind timing process to its selected core")
    record = {"slurm_allowed_cpus": ",".join(map(str, allowed)), "timing_cpu": cpu, "selection": method,
              "gpu_logical_ordinal": 0, "gpu_pci": pci, "gpu_numa_node": node,
              "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", "")}
    with (directory / "binding.tsv").open("x", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=record, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(record)


def prepare(directory, selected, cli):
    published = {}
    with (ROOT / "results/repeated/divergence-gh200/divergence.csv").open(newline="") as source:
        for row in csv.DictReader(source):
            published[row["path"]] = (row["sha256"], int(row["distance"]))
    conditions = {}
    for group in selected:
        length, rate, reps = GROUPS[group]
        pair_groups = conditions.setdefault((length, rate), {})
        pair_groups.update((rep, group) for rep in reps)
    manifest = []
    for (length, rate), pair_groups in conditions.items():
        generate(directory, [length], [Decimal(rate)], max(pair_groups), MASTER_SEED)
        with (directory / "manifest.tsv").open(newline="") as source:
            for row in csv.DictReader(source, delimiter="\t"):
                if int(row["rep"]) in pair_groups:
                    manifest.append((pair_groups[int(row["rep"])], row))
    inputs = []
    for group, row in manifest:
        path = directory / row["path"]
        digest = sha256(path)
        result = subprocess.check_output(
            [cli, "pair", "--algorithm", "edlib", "--threads", "1", "--input", path], text=True)
        score = int(json.loads(result)["distance"])
        if group != "held-out" and published.get(row["path"]) != (digest, score):
            raise ValueError(f"Historical FASTA or score mismatch: {path.name}")
        inputs.append({"group": group, **row, "sha256": digest, "distance": score})
    if len(inputs) != 5 * len(selected):
        raise ValueError("Missing input in a requested group")
    return inputs


def run(binary, method, item, directory):
    algorithm = "nvpd-cuda" if method in ("original", "candidate") else method
    output = subprocess.check_output(
        [binary, directory / item["path"], algorithm, str(item["distance"])],
        input="run\nrun\n", text=True).splitlines()
    if len(output) != 3 or output[0] != "ready":
        raise ValueError(f"Unexpected repeat_kernel output for {method}: {output}")
    times = []
    for line in output[1:]:
        m, n, score, seconds, rss, device = line.split(",")
        seconds = float(seconds)
        if ((int(m), int(n), int(score)) != (int(item["length"]), int(item["mutant_length"]), item["distance"])
                or not math.isfinite(seconds) or seconds <= 0):
            raise ValueError(f"Incorrect score, dimensions or time for {method}: {line}")
        times.append(seconds)
    return {"distance": int(score), "warmup_seconds": times[0], "seconds": times[1],
            "rss_kib": int(rss), "device_bytes": int(device)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("original", "candidate", "editdist", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--blocks", type=int, default=8, help="positive multiple of four")
    parser.add_argument("--groups", nargs="+", choices=GROUPS, default=DEFAULT_GROUPS)
    args = parser.parse_args()
    if args.blocks < 4 or args.blocks % 4 or len(set(args.groups)) != len(args.groups):
        parser.error("Supply a multiple of four blocks and unique group names")
    binaries = {m: (args.original if m == "original" else args.candidate).resolve() for m in METHODS}
    cli = args.editdist.resolve()
    os.environ.update(OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_PLACES="cores",
                      OMP_PROC_BIND="close", OMP_WAIT_POLICY="ACTIVE")
    args.output.mkdir(parents=True, exist_ok=True)
    pin_timing_process(args.output)
    with (args.output / "binaries.tsv").open("x", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["job", "program", "sha256"])
        for name, path in (("original", binaries["original"]), ("candidate", binaries["candidate"]), ("edlib-oracle", cli)):
            writer.writerow([os.environ.get("SLURM_JOB_ID", "local"), name, sha256(path)])
    with TemporaryDirectory(prefix="nvpd-focus-") as temp:
        directory = Path(temp)
        inputs = prepare(directory, args.groups, cli)
        with (args.output / "inputs.tsv").open("x", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=inputs[0], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(inputs)
        rng = random.Random(67431)
        with (args.output / "trials.csv").open("x", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=TRIAL_FIELDS, lineterminator="\n")
            writer.writeheader()
            for block in range(1, args.blocks + 1):
                order = list(enumerate(inputs))
                rng.shuffle(order)
                for i, item in order:
                    # Four Williams orders balance positions and consecutive method pairs.
                    for position, offset in enumerate((0, 1, 3, 2), start=1):
                        method = METHODS[(i + block - 1 + offset) % 4]
                        result = run(binaries[method], method, item, directory)
                        writer.writerow({"group": item["group"], "rep": item["rep"], "block": block,
                                         "order": position, "method": method, **result})
                    target.flush()
                print(f"Finished block {block}/{args.blocks}", flush=True)


if __name__ == "__main__":
    main()
