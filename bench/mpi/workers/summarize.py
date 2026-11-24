"""Block-paired Kokkos/OpenMP times and strong scaling on one and two nodes."""
import csv
import math
import statistics
import sys
from pathlib import Path


def summarize(folder):
    with (folder / "trials.csv").open() as f:
        rows = list(csv.DictReader(f))
    observations = {(int(r["block"]), r["schedule"], int(r["nodes"]), r["backend"]):r for r in rows}
    if len(rows) != 100 or len(observations) != 100 or len({r["cells"] for r in rows}) != 1:
        raise ValueError("expected one complete fixed-work experiment")
    with (folder / "summary.csv").open("w") as out:
        writer = csv.writer(out)
        writer.writerow(("schedule","nodes","backend","metric","n","mean","sd","low95","high95"))
        for schedule in sorted({r["schedule"] for r in rows}):
            for nodes in (1,2):
                for backend in ("openmp","kokkos"):
                    samples = [observations[b,schedule,nodes,backend] for b in range(1,6)]
                    seconds = [float(r["seconds"]) for r in samples]
                    metrics = {"seconds":seconds, "peak_rank_rss_kib":[int(r["peak_rank_rss_kib"]) for r in samples]}
                    if backend == "kokkos":
                        metrics["speedup_over_openmp"] = [float(observations[b,schedule,nodes,"openmp"]["seconds"])/t
                                                         for b,t in enumerate(seconds,1)]
                    if nodes == 2:
                        metrics["strong_efficiency"] = [float(observations[b,schedule,1,backend]["seconds"])/(2*t)
                                                       for b,t in enumerate(seconds,1)]
                    for metric,values in metrics.items():
                        mean,sd = statistics.mean(values),statistics.stdev(values)
                        half = 2.7764451052*sd/math.sqrt(5)
                        writer.writerow((schedule,nodes,backend,metric,5,mean,sd,mean-half,mean+half))


if __name__ == "__main__":
    summarize(Path(sys.argv[1]))
