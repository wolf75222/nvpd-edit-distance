#!/usr/bin/env python3
"""Validate the 1M completion campaign against both libraries and the original measurements."""
import argparse
import csv
import math
from pathlib import Path
from statistics import mean, median, stdev

RATES = (0.02, 0.05, 0.1)
DENSE = ("diagonal", "nvpd-omp")
LIBS = ("edlib", "wfa2")
REVISION = "f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83"


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows, sep=","):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=rows[0], delimiter=sep, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def index(rows):
    result = {}
    for row in rows:
        key = row["algorithm"], float(row["rate"]), int(row["rep"])
        if key in result:
            raise ValueError(f"Duplicate trial at {key}")
        if (row["status"] != "ok" or int(row["length"]) != 1000000 or
                int(row["m"]) != 1000000 or int(row["trial"]) != 1 or
                int(row["threads"]) != (8 if key[0] in DENSE else 1) or
                not math.isfinite(float(row["seconds"])) or float(row["seconds"]) <= 0):
            raise ValueError(f"Unexpected trial parameters at {key}")
        result[key] = row
    return result


def check_pair(x, y):
    for field in ("input_sha256", "path", "seed", "m", "n", "distance"):
        if x[field] != y[field]:
            raise ValueError(f"{field} mismatch at {x['path']}")


def cpus(text):
    ids = set()
    for part in text.split(","):
        bounds = list(map(int, part.split("-")))
        ids.update(range(bounds[0], bounds[-1] + 1))
    return ids


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job", type=Path, help="directory containing raw/*.csv and rank CPU masks")
    parser.add_argument("historical", type=Path, help="original results/romeo directory")
    args = parser.parse_args()
    cells = {(alg, rate, rep) for alg in DENSE + LIBS for rate in RATES for rep in range(1, 6)}
    raw = args.job / "raw"
    timed, placement = [], []
    seen_cpus = {}
    for path in sorted(raw.glob("rep-*.csv")):
        data = read(path)
        if len(data) != 1:
            raise ValueError(f"Expected a single trial in {path}")
        row = data[0]
        mask = path.with_suffix(".cpus").read_text().strip()
        ids = cpus(mask)
        if len(ids) != 8:
            raise ValueError(f"Expected eight bound CPU IDs at {path}, found {mask}")
        key = row["rep"], row["host"]
        used = seen_cpus.setdefault(key, set())
        if used.intersection(ids):
            raise ValueError(f"Overlapping CPU masks in wave {key}")
        used.update(ids)
        placement.append(dict(job=row["job"], wave=row["rep"], host=row["host"],
                              algorithm=row["algorithm"], rate=row["rate"], cpu_ids=mask))
        timed.append(row)
    if len(seen_cpus) != 10 or any(len(group) != 24 for group in seen_cpus.values()):
        raise ValueError("Expected two groups of 24 disjoint cores for each of five waves")
    refs = [r for p in sorted(raw.glob("refs-*.csv")) for r in read(p)]
    new = index(timed + refs)
    if set(new) != cells:
        raise ValueError(f"Incomplete campaign: missing {cells.difference(new)}, extra {set(new).difference(cells)}")

    old_lib = index([r for r in read(args.historical / "cpu-837264.csv")
                     if r["length"] == "1000000" and r["algorithm"] in LIBS])
    original_dp = [r for j in (837266, 837267, 837268)
                   for r in read(args.historical / f"cpu-{j}.csv")]
    old_dense = index([r for r in original_dp if r["algorithm"] in DENSE])
    old_seq = index([r for r in original_dp if r["algorithm"] == "seq"])
    if len(old_lib) != 30 or len(old_dense) != 6 or len(old_seq) != 3:
        raise ValueError("Original library or control rows missing")
    original_hashes = {r["binary_sha256"] for table in (old_lib, old_dense, old_seq) for r in table.values()}
    run_hashes = {r["binary_sha256"] for r in new.values()}
    if len(run_hashes) != 1:
        raise ValueError("More than one binary was used in this run")
    same_binary = run_hashes == original_hashes

    if len({r["input_sha256"] for r in new.values()}) != 15:
        raise ValueError("Expected 15 distinct pair FASTAs")
    cover, summary = [], []
    for alg in DENSE:
        for rate in RATES:
            rows = []
            for rep in range(1, 6):
                row = new[alg, rate, rep]
                for lib in LIBS:
                    check_pair(row, old_lib[lib, rate, rep])
                    check_pair(row, new[lib, rate, rep])
                old = old_dense.get((alg, rate, rep))
                serial = old_seq.get(("seq", rate, rep))
                if old:
                    check_pair(row, old)
                if serial:
                    check_pair(row, serial)
                t = float(row["seconds"])
                rows.append(row)
                cover.append(dict(algorithm=alg, rate=rate, rep=rep, input_sha256=row["input_sha256"],
                                  distance=row["distance"], m=row["m"], n=row["n"],
                                  binary_matches_original=int(same_binary),
                                  original_job=old["job"] if old else "",
                                  original_seconds=old["seconds"] if old else "", current_job=row["job"],
                                  current_seconds=row["seconds"], host=row["host"], threads=row["threads"],
                                  current_over_original=t / float(old["seconds"]) if old else "",
                                  historical_seq_seconds=serial["seconds"] if serial else "",
                                  current_speedup_vs_historical_seq=float(serial["seconds"])/t if serial else "",
                                  current_efficiency_vs_historical_seq=float(serial["seconds"])/t/8 if serial else "",
                                  coverage="repeated_control" if old else "previously_missing"))
            ts = [float(x["seconds"]) for x in rows]
            gcups = [int(x["m"])*int(x["n"])/float(x["seconds"])/1e9 for x in rows]
            summary.append(dict(algorithm=alg, rate=rate, revision=REVISION, threads=8, pairs=5,
                                min_seconds=min(ts), median_seconds=median(ts), mean_seconds=mean(ts),
                                sample_sd_seconds=stdev(ts), max_seconds=max(ts), median_gcups=median(gcups),
                                max_peak_rss_kib=max(int(x["peak_rss_kib"]) for x in rows),
                                rep1_current_over_original=ts[0]/float(old_dense[alg, rate, 1]["seconds"])))
    key = lambda r: (r["algorithm"], float(r["rate"]), int(r["rep"]))
    write(args.job / "trials.csv", sorted(timed, key=key))
    write(args.job / "references.csv", sorted(refs, key=key))
    write(args.job / "coverage.tsv", cover, "\t")
    write(args.job / "summary.csv", summary)
    write(args.job / "placement.tsv", placement, "\t")
    print("Validated 30 eight-thread timings, 30 fresh library calls, 30 original library results, 6 repeated controls and 3 original sequential references.")


if __name__ == "__main__":
    main()
