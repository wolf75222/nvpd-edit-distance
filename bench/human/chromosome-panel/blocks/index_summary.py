"""Validate both lanes and aggregate the nine independent, paired process blocks per case."""
import argparse
import csv
import hashlib
from math import isfinite, sqrt
from pathlib import Path
from statistics import mean, stdev

from index import protocol


METRICS = ("seconds", "peak_rss_kib", "allocations", "requested_bytes", "instructions", "cycles",
           "L1_loads", "L1_misses")
T95_DF8 = 2.306004135204166


def statistics(values):
    if len(values) != 9 or not all(isfinite(x) and x > 0 for x in values):
        raise ValueError("expected nine positive observations")
    mu, sd = mean(values), stdev(values)
    half = T95_DF8 * sd / sqrt(len(values))
    return dict(n=len(values), mean=mu, sd=sd, low95=mu-half, high95=mu+half)


def csv_rows(path, delimiter=","):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def save(path, rows, delimiter=","):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--htslib-indexes", type=Path, help="also compare all six files byte for byte with samtools")
    args = p.parse_args()
    raw = [row for f in sorted(args.raw.glob("index-*.csv")) for row in csv_rows(f)]
    if len(raw) != 108 or len({r["job"] for r in raw}) != 1:
        raise ValueError("expected 108 timed processes in one allocation")
    if len({(r["lane"], r["cpu"]) for r in raw}) != 2 or len({r["cpu"] for r in raw}) != 2:
        raise ValueError("expected two disjoint fixed CPU placements")
    summary, paired, entries = [], [], []
    for case in protocol.READ_CASES:
        rows = [r for r in raw if r["case"] == case]
        keyed = {(r["method"], int(r["block"])): r for r in rows}
        if len(rows) != 18 or len(keyed) != 18 or len({(r["lane"], r["cpu"]) for r in rows}) != 1:
            raise ValueError(f"invalid blocks or affinity for {case}")
        for block in range(1, 10):
            if {keyed[m, block]["order"] for m in ("blocks", "indexed")} != {"1", "2"}:
                raise ValueError("incomplete paired block")
        if any(float(r["pmu_min_running_pct"]) < 99 for r in rows):
            raise ValueError("multiplexed PMU")
        index = args.raw / (case + ".fa.fai")
        if args.htslib_indexes and index.read_bytes() != (args.htslib_indexes / index.name).read_bytes():
            raise ValueError(f"samtools index mismatch for {case}")
        for field in ("input_sha256", "normalized_sha256", "fai_sha256"):
            if len({r[field] for r in rows}) != 1:
                raise ValueError("input or decoded identity changed")
        if hashlib.sha256(index.read_bytes()).hexdigest() != rows[0]["fai_sha256"]:
            raise ValueError("collected FAI digest differs")
        with index.open() as f:
            for name, length, offset, bases, width in csv.reader(f, delimiter="\t"):
                entries.append(dict(case=case, name=name, length=length, offset=offset,
                                    linebases=bases, linewidth=width))
        for metric in METRICS:
            series = {}
            for method in ("blocks", "indexed"):
                values = [float(keyed[method, i][metric]) for i in range(1, 10)]
                series[method] = values
                summary.append(dict(case=case, method=method, metric=metric, **statistics(values)))
            ratios = [a/b for a, b in zip(series["blocks"], series["indexed"])]
            paired.append(dict(case=case, numerator="blocks", denominator="indexed", metric=metric,
                               **statistics(ratios)))
    args.out.mkdir(parents=True, exist_ok=True)
    save(args.out/"index-trials.csv", raw)
    save(args.out/"index-summary.csv", summary)
    save(args.out/"index-paired.csv", paired)
    save(args.out/"index-entries.tsv", entries, "\t")


if __name__ == "__main__":
    main()
