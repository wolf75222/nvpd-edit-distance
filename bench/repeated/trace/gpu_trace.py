"""Extract the Nsight per-operation rows and summarize the seven independent processes."""
import csv
from pathlib import Path
import statistics
import sys


src, out = map(Path, sys.argv[1:])
reports = [("gpu_kern", "kernel", "Instances", "Name"),
           ("gpu_mem_time", "copy", "Count", "Operation"),
           ("api", "api", "Num Calls", "Name")]
ops, captures = [], []
for block in range(1, 8):
    order = ("baseline", "candidate") if block % 2 else ("candidate", "baseline")
    for position, variant in enumerate(order, 1):
        prefix = str(src / f"{block}-{variant}")
        lines = [r for r in csv.reader(Path(prefix + ".log").read_text().splitlines()) if r and r[0] == variant]
        if len(lines) != 2 or [r[4] for r in lines] != ["cold", "warm"]:
            raise ValueError(f"Missing warmup or captured invocation in {prefix}")
        for r in lines:
            if tuple(int(r[i]) for i in (2, 3, 5, 7, 9, 10)) != (5000, 4992, block, 0, 99, 120024):
                raise ValueError(f"Unexpected dimensions, run, trace mode, score or device bytes in {prefix}")
        row = dict(block=block, position=position, variant=variant, m=5000, n=4992, distance=99, device_bytes=int(lines[-1][10]),
                   wall_ms=1000 * float(lines[-1][8]))
        for report, domain, count_key, name_key in reports:
            with Path(f"{prefix}_cuda_{report}_sum.csv").open(newline="") as stream:
                data = list(csv.DictReader(stream))
            if not data:
                raise ValueError(f"Empty {domain} domain in {prefix}")
            for r in data:
                ops.append(dict(block=block, variant=variant, domain=domain, name=r[name_key],
                                count=int(r[count_key]), total_ns=int(r["Total Time (ns)"])))
            row[domain + "_ms"] = sum(int(r["Total Time (ns)"]) for r in data) / 1e6
            if domain == "kernel":
                row["launches"] = sum(int(r[count_key]) for r in data)
        if row["launches"] != (4993 if variant == "baseline" else 625):
            raise ValueError(f"Unexpected number of kernels: {row}")
        captures.append(row)


def write(name, rows):
    with (out / (name + ".csv")).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


write("operations", ops)
write("captures", captures)
rows = []
for variant in ("baseline", "candidate"):
    for metric in ("wall_ms", "kernel_ms", "copy_ms", "api_ms", "launches"):
        sample = [row[metric] for row in captures if row["variant"] == variant]
        mean, sd = statistics.mean(sample), statistics.stdev(sample)
        half = 2.446911851 * sd / len(sample)**0.5  # t(0.975, 6)
        rows.append(dict(variant=variant, metric=metric, processes=len(sample), mean=mean, sd=sd,
                         low95=mean-half, high95=mean+half))
write("summary", rows)
