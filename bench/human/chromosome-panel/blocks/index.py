"""Build standard five-field FAI files and compare two modes of one Linux binary."""
import csv
import hashlib
import os
from pathlib import Path
import random
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import measure as protocol


def make_index(fasta, index):
    """Write NAME, LENGTH, OFFSET, LINEBASES, LINEWIDTH for fixed-width panel FASTA."""
    name = None
    offset = total = length = bases = width = previous = 0
    with fasta.open("rb") as src, index.open("x") as dst:
        for line in src:
            total += len(line)
            if line.startswith(b">"):
                if name is not None:
                    print(name, length, offset, bases, width, sep="\t", file=dst)
                name = line[1:].split()[0].decode("ascii")
                length, bases, width, previous, offset = 0, 0, 0, 0, total
            else:
                n = len(line.rstrip(b"\r\n"))
                if name is None or not n:
                    raise ValueError("expected nonempty reference sequence lines")
                if not length:
                    bases, width = n, len(line)
                elif previous != width or n > bases:
                    raise ValueError("nonterminal line changes width")
                length += n
                previous = len(line)
        if name is None:
            raise ValueError("empty panel FASTA")
        print(name, length, offset, bases, width, sep="\t", file=dst)


def sequence_digest(binary, fasta, extra):
    with subprocess.Popen([binary, fasta, "check", *extra], stdout=subprocess.PIPE) as proc:
        proc.stdout.readline()  # length, score placeholder, time and RSS
        digest = hashlib.sha256()
        for chunk in iter(lambda: proc.stdout.read(1024 * 1024), b""):
            digest.update(chunk)
        if proc.wait():
            raise RuntimeError("failed to normalize FASTA")
        return digest.hexdigest()


def counters(binary, path, extra, case):
    lines = subprocess.check_output([binary, path, "profile", *extra], text=True, timeout=240).splitlines()
    if len(lines) != 8:
        raise ValueError("missing allocation or PMU measurement")
    for line in lines[:2]:
        m, n = line.split(",")[:2]
        if (m, n) != (case["m"], case["n"]):
            raise ValueError("profile length mismatch")
    values = {}
    for line in lines[2:]:
        name, value, *pct = line.split(",")
        if pct:
            if float(pct[0]) < 99:
                raise ValueError("PMU group was multiplexed")
            values["pmu_min_running_pct"] = min(float(pct[0]), values.get("pmu_min_running_pct", 100))
        values[name] = int(value)
    values["L1_loads"] = values.pop("L1-dcache-loads")
    values["L1_misses"] = values.pop("L1-dcache-load-misses")
    return values


def main():
    data, binary, out = map(Path, sys.argv[1:])
    lane = int(os.environ["SLURM_PROCID"])
    lanes = int(os.environ["SLURM_NTASKS"])
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) != 1:
        raise ValueError(f"expected one physical CPU, got {cpus}")
    with (data / "manifest.tsv").open() as stream:
        manifest = {row["case"]: row for row in csv.DictReader(stream, delimiter="\t")}
    common = dict(job=os.environ["SLURM_JOB_ID"], lane=lane, cpu=cpus[0])
    fields = ("job lane cpu case input_sha256 normalized_sha256 fai_sha256 method block order "
              "warm_seconds seconds peak_rss_kib instructions cycles L1_loads L1_misses "
              "pmu_min_running_pct allocations requested_bytes")
    rng = random.Random(83291 + lane)
    with (out / f"index-{lane}.csv").open("x") as stream:
        sink = protocol.writer(stream, fields)
        for name in protocol.READ_CASES[lane::lanes]:
            fasta = data / (name + ".fa")
            index = out / (name + ".fa.fai")
            make_index(fasta, index)
            with fasta.open("rb") as src:
                file_sha = hashlib.sha256()
                for chunk in iter(lambda: src.read(1024 * 1024), b""):
                    file_sha.update(chunk)
            if file_sha.hexdigest() != manifest[name]["input_sha256"]:
                raise ValueError("input SHA-256 mismatch")
            expected = protocol.canonical_file_digest(fasta)
            modes = {"blocks": [], "indexed": [index]}
            for extra in modes.values():
                if sequence_digest(binary, fasta, extra) != expected:
                    raise ValueError("normalized FASTA differs from independent SHA-256")
            identity = dict(common, case=name, input_sha256=file_sha.hexdigest(), normalized_sha256=expected,
                            fai_sha256=hashlib.sha256(index.read_bytes()).hexdigest())
            cases = list(modes)
            for block in range(1, 10):
                rng.shuffle(cases)
                for order, mode in enumerate(cases, 1):
                    args = modes[mode]
                    timing = protocol.run_trial([binary, fasta, "reader", *args], manifest[name])
                    del timing["distance"]
                    events = counters(binary, fasta, args, manifest[name])
                    sink.writerow(dict(identity, method=mode, block=block, order=order, **timing, **events))
                    stream.flush()
            print(f"{name}: nine blocks, complete FASTA digest verified", flush=True)


if __name__ == "__main__":
    main()
