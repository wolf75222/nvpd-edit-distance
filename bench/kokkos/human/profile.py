"""Summarize physical CUDA launches in one separately profiled, warmed backend call."""

import argparse
import csv
from pathlib import Path
import sqlite3
import subprocess

from measure import write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("binary", "inputs", "out", "temp"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    refs = {r["case"]: r for r in csv.DictReader((a.out / "references.csv").open())}
    rows = []
    totals = []
    for case in ("collinear-16384", "insertion-16384", "collinear-65536"):
        for alg in ("native", "kokkos-rows", "kokkos-opt"):
            stem = a.temp / f"{case}-{alg}"
            subprocess.run(
                [
                    "nsys",
                    "profile",
                    "--trace=cuda",
                    "--sample=none",
                    "--cpuctxsw=none",
                    "--capture-range=cudaProfilerApi",
                    "--capture-range-end=stop",
                    "--output",
                    stem,
                    a.binary,
                    alg,
                    a.inputs / (case + ".fa"),
                    "1",
                    refs[case]["distance"],
                ],
                check=True,
            )
            db = stem.with_suffix(".sqlite")
            subprocess.run(
                ["nsys", "export", "--type=sqlite", "-o", db, stem.with_suffix(".nsys-rep")], check=True
            )
            with sqlite3.connect(db) as conn:
                sql = (
                    "SELECT s.value, count(*), sum(k.end-k.start) FROM CUPTI_ACTIVITY_KIND_KERNEL k "
                    "JOIN StringIds s ON s.id=k.demangledName GROUP BY s.value ORDER BY count(*) DESC"
                )
                kernels = conn.execute(sql).fetchall()
                transfers = conn.execute(
                    "SELECT count(*), sum(bytes),sum(end-start) FROM CUPTI_ACTIVITY_KIND_MEMCPY"
                ).fetchone()
            if not kernels:
                raise RuntimeError("empty GPU trace")
            common = dict(case=case, algorithm=alg, input_sha256=refs[case]["input_sha256"])
            for index, (name, count, ns) in enumerate(kernels):
                rows.append(
                    dict(**common, group=index + 1, kernel=name, launches=count, gpu_seconds=ns / 1e9)
                )
            totals.append(
                dict(
                    **common,
                    kernel_launches=sum(x[1] for x in kernels),
                    kernel_seconds=sum(x[2] for x in kernels) / 1e9,
                    copies=transfers[0],
                    copy_bytes=transfers[1],
                    copy_seconds=transfers[2] / 1e9,
                )
            )
            write(a.out / "kernel_groups.csv", rows)
            write(a.out / "profile.csv", totals)
    print("Nine capture ranges completed and their edit distances matched both libraries.")


if __name__ == "__main__":
    main()
