"""Collect CUDA activity for one complete NvPD score on an identified FASTA pair."""
import argparse
import csv
from pathlib import Path
import sqlite3
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("score", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--algorithm", default="nvpd-cuda", choices=("nvpd-cuda", "nvpd-cuda-persistent"))
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        stem = Path(directory) / "activity"
        subprocess.run(["nsys", "profile", "--sample=none", "--cpuctxsw=none", "--trace=cuda",
                        "-o", stem, args.binary.resolve(), args.fasta, args.algorithm, str(args.score)],
                       input="run\n", text=True, check=True)
        db = stem.with_suffix(".sqlite")
        subprocess.run(["nsys", "export", "-t", "sqlite", "-o", db, stem.with_suffix(".nsys-rep")],
                       check=True)
        with sqlite3.connect(db) as conn, args.output.open("x", newline="") as out:
            writer = csv.writer(out, lineterminator="\n")
            writer.writerow(("activity", "name", "count", "total_seconds", "total_bytes"))
            writer.writerows(conn.execute("""
                SELECT 'kernel', names.value, COUNT(*), SUM(k.end-k.start)/1e9, NULL
                FROM CUPTI_ACTIVITY_KIND_KERNEL k JOIN StringIds names ON names.id=k.shortName
                GROUP BY names.value ORDER BY SUM(k.end-k.start) DESC"""))
            writer.writerows(conn.execute("""
                SELECT 'copy', copyKind, COUNT(*), SUM(end-start)/1e9, SUM(bytes)
                FROM CUPTI_ACTIVITY_KIND_MEMCPY GROUP BY copyKind ORDER BY copyKind"""))


if __name__ == "__main__":
    main()
