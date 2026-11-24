"""Five matched process blocks per frozen human chromosome interval."""

import argparse
import csv
import hashlib
import os
from pathlib import Path
import random
import sqlite3
import subprocess


def profiler(binary, inputs, refs, out, tmp, local, classes):
    methods = {"native": "next_row", "fused": "pair_rows"}
    methods.update({"local": "local_row_tile"} if local else {"shuffle": "pair_rows"})
    if classes:
        methods = {"local": "local_row_tile", "classes": "class_row_tile"}
    with (out / "activity.csv").open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            (
                "case",
                "method",
                "kernel_launches",
                "kernel_seconds",
                "copies",
                "copy_bytes",
                "copy_seconds",
            )
        )
        for case in ("collinear-16384", "insertion-524288" if classes else "insertion-131072"):
            for method in methods:
                stem = tmp / f"{case}-{method}"
                cmd = [binary, method, inputs / (case + ".fa"), refs[case]["distance"], "profile"]
                subprocess.run(
                    [
                        "nsys",
                        "profile",
                        "--trace=cuda",
                        "--sample=none",
                        "--cpuctxsw=none",
                        "--capture-range=cudaProfilerApi",
                        "--capture-range-end=stop",
                        "-o",
                        stem,
                        *cmd,
                    ],
                    stdout=subprocess.DEVNULL,
                    check=True,
                )
                db = stem.with_suffix(".sqlite")
                subprocess.run(
                    ["nsys", "export", "--type=sqlite", "-o", db, stem.with_suffix(".nsys-rep")],
                    stdout=subprocess.DEVNULL,
                    check=True,
                )
                with sqlite3.connect(db) as conn:
                    nk, ns = conn.execute(
                        "SELECT count(*),sum(end-start) FROM CUPTI_ACTIVITY_KIND_KERNEL"
                    ).fetchone()
                    nc, nb, ncs = conn.execute(
                        "SELECT count(*),sum(bytes),sum(end-start) FROM CUPTI_ACTIVITY_KIND_MEMCPY"
                    ).fetchone()
                writer.writerow((case, method, nk, ns / 1e9, nc, nb, ncs / 1e9))
                for path in (db, stem.with_suffix(".nsys-rep")):
                    path.unlink()
    for method, kernel in methods.items():
        case = "collinear-16384" if classes else "insertion-131072"
        log = out / f"ncu-{method}.txt"
        cmd = [binary, method, inputs / (case + ".fa"), refs[case]["distance"], "profile"]
        result = subprocess.run(
            [
                "ncu",
                "--profile-from-start",
                "off",
                "--kernel-name",
                f"regex:{kernel}",
                "--launch-count",
                "1",
                "--metrics",
                "dram__bytes_read.sum,dram__bytes_write.sum,smsp__inst_executed.sum,gpu__time_duration.sum",
                "--csv",
                "--page",
                "raw",
                "--log-file",
                log,
                *cmd,
            ],
            capture_output=True,
        )
        if result.returncode and "ERR_NVGPUCTRPERM" not in log.read_text():
            result.check_returncode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "inputs", "out", "tmp", "references"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--local", action="store_true")
    parser.add_argument("--classes", action="store_true")
    a = parser.parse_args()
    with a.references.open() as f:
        refs = {r["case"]: r for r in csv.DictReader(f)}
    chosen = [r for r in refs.values() if int(r["m"]) in (16384, 65536, 131072, 524288)]
    rng = random.Random(9126)
    with (a.out / "trials.csv").open("x") as stream:
        w = csv.writer(stream, lineterminator="\n")
        w.writerow(
            (
                "job",
                "case",
                "input_sha256",
                "block",
                "order",
                "method",
                "m",
                "n",
                "distance",
                "seconds",
                "peak_rss_kib",
                "device_bytes",
            )
        )
        for ref in chosen:
            case = ref["case"]
            path = a.inputs / (case + ".fa")
            if hashlib.sha256(path.read_bytes()).hexdigest() != ref["input_sha256"]:
                raise ValueError("FASTA digest mismatch")
            methods = ["native", "fused", "local" if a.local else "shuffle", "edlib", "wfa2"]
            if a.classes:
                methods = ["native", "fused", "local", "classes", "edlib"]
            for block in range(1, 6):
                rng.shuffle(methods)
                for order, method in enumerate(methods, 1):
                    line = (
                        subprocess.check_output(
                            [a.binary, method, path, ref["distance"], "time"], text=True
                        )
                        .strip()
                        .split(",")
                    )
                    if line[:3] != [ref[k] for k in ("m", "n", "distance")] or float(line[3]) <= 0:
                        raise ValueError("incorrect score, dimensions or timer")
                    w.writerow(
                        (
                            os.environ["SLURM_JOB_ID"],
                            case,
                            ref["input_sha256"],
                            block,
                            order,
                            method,
                            *line,
                        )
                    )
                    stream.flush()
            print(case, "complete", flush=True)
    profiler(a.binary, a.inputs, refs, a.out, a.tmp, a.local, a.classes)


if __name__ == "__main__":
    main()
