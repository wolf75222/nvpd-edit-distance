"""Five randomized complete blocks on one binding, each with a fresh timed process."""
import csv
import os
from pathlib import Path
import random
import subprocess
import sys


def main():
    executable, dataset, output = sys.argv[1:]
    root = Path(dataset)
    with (root / "manifest.tsv").open() as f:
        sizes = {(row["chrom"], row["version"]): int(row["bases"])
                 for row in csv.DictReader(f, delimiter="\t")}
    chroms = sorted({chrom for chrom, _ in sizes})
    matrix = [(c, m, t) for c in chroms for m in ("split", "sparse") for t in (1, 2, 4, 8)]
    matrix.extend((c, m, 1) for c in chroms for m in ("edlib", "wfa2", "biwfa"))
    rng = random.Random(261125)
    scores = {}
    with open(output, "x", buffering=1) as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("job", "block", "chrom", "method", "threads", "m", "n", "score",
                         "warm_s", "seconds", "rss_kib"))
        for block in range(1, 6):
            rng.shuffle(matrix)
            for chrom, method, threads in matrix:
                args = [executable, root / f"{chrom}.v1.1.fa", root / f"{chrom}.v1.2.fa", method]
                env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PLACES="cores",
                           OMP_PROC_BIND="close", OMP_DYNAMIC="false")
                result = subprocess.check_output(args, text=True, env=env, timeout=120)
                warm, timed = (row.split(",") for row in result.splitlines())
                if warm[:3] != timed[:3] or tuple(map(int,timed[:2])) != (sizes[chrom,"1.1"],sizes[chrom,"1.2"]):
                    raise ValueError(f"{chrom}: timing records disagree")
                if int(timed[2]) != scores.setdefault(chrom, int(timed[2])):
                    raise ValueError(f"{chrom}: scores disagree")
                writer.writerow((env["SLURM_JOB_ID"], block, chrom, method, threads,
                                 *timed[:3], warm[3], *timed[3:]))


    with Path(output).with_name("counters.csv").open("x", buffering=1) as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("job","block","chrom","method","event","count","running_pct"))
        methods = [(c, m) for c in chroms for m in ("sparse", "split")]
        env = dict(os.environ, OMP_NUM_THREADS="1", OMP_PLACES="cores", OMP_PROC_BIND="close")
        for block in range(1, 6):
            rng.shuffle(methods)
            for chrom, method in methods:
                args = [executable,root / f"{chrom}.v1.1.fa",root / f"{chrom}.v1.2.fa",method,"counters"]
                lines = subprocess.check_output(args, text=True, env=env, timeout=120).splitlines()
                if int(lines[1].split(",")[2]) != scores[chrom] or len(lines) != 6:
                    raise ValueError("counter trial disagrees with score benchmark")
                for raw in lines[2:]:
                    event, count, pct = raw.split(",")
                    if float(pct) < 99:
                        raise ValueError(f"multiplexed hardware event {event}")
                    writer.writerow((env["SLURM_JOB_ID"],block,chrom,method,event,count,pct))


if __name__ == "__main__":
    main()
