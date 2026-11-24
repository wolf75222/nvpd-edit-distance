#!/usr/bin/env python3
"""Check every block and combine Slurm lanes for the CHM13 panel."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev


def read(path, delimiter=","):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(values):
    if len(values) != 5:
        raise ValueError("expected five measurements")
    mu, sd = mean(values), stdev(values)
    half = 2.7764451051977987 * sd / 5**.5
    return dict(n=5, mean=mu, sd=sd, low95=mu-half, high95=mu+half)


def collect_blocks(lanes, out):
    from measure import READ_CASES

    records = [item for path in sorted(lanes.glob("blocks-*.csv")) for item in read(path)]
    data = {(r["case"], r["method"], int(r["block"])): r for r in records}
    methods = ("portable", "blocks")
    expected = [(case, method, block) for case in READ_CASES for block in range(1, 6) for method in methods]
    if len(records) != len(expected) or data.keys() != set(expected):
        raise ValueError("incomplete or duplicate SIMD block experiment")
    for case in READ_CASES:
        if len({data[case, method, block]["cpu"] for method in methods for block in range(1, 6)}) != 1:
            raise ValueError(f"worker migrated between trials on {case}")
    out.mkdir(exist_ok=True, parents=True)
    write(out / "blocks.csv", [data[key] for key in expected])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("lanes", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--reader-blocks", action="store_true", help="combine the explicit SIMD comparison")
    opt = p.parse_args()
    if opt.reader_blocks:
        collect_blocks(opt.lanes, opt.out)
        return
    manifest = {r["case"]: r for r in read(opt.lanes / "manifest.tsv", "\t")}
    if len(manifest) != 18:
        raise ValueError("expected 18 complete chromosomes")
    opt.out.mkdir(parents=True, exist_ok=True)
    (opt.out / "manifest.tsv").write_bytes((opt.lanes / "manifest.tsv").read_bytes())
    data = {}
    for kind, size in (("trials", 180), ("reader", 90), ("counters", 360), ("checks", 18)):
        files = sorted(opt.lanes.glob(f"{kind}-*.csv"))
        rows = [row for f in files for row in read(f)]
        keys = [(r["case"], r["method"], r.get("block"), r.get("event")) for r in rows]
        if len(files) != 8 or len(keys) != len(set(keys)) or len(rows) != size:
            raise ValueError(f"incomplete or duplicate {kind}")
        data[kind] = sorted(rows, key=lambda r: (list(manifest).index(r["case"]), r["method"],
                                                  int(r.get("block", 0)), r.get("event", "")))
        write(opt.out / f"{kind}.csv", data[kind])
    checksets = defaultdict(set)
    for row in data["checks"]:
        checksets[row["case"]].add(row["normalized_sha256"])
    if any(len(s) != 1 for s in checksets.values()) or len(checksets) != 6:
        raise ValueError("FASTA readers differ")
    scores = defaultdict(set)
    summaries, ratios = [], []
    for kind in ("trials", "reader"):
        groups = defaultdict(dict)
        for row in data[kind]:
            case, method, block = row["case"], row["method"], int(row["block"])
            if kind == "trials":
                scores[case].add(int(row["distance"]))
            groups[case, method][block] = row
        methods = ("edlib", "wfa2") if kind == "trials" else ("initial", "scalar", "simd")
        for case in {key[0] for key in groups}:
            for method in methods:
                items = groups[case, method]
                if set(items) != set(range(1,6)): raise ValueError((case,method,"blocks"))
                for measure in ("seconds", "peak_rss_kib"):
                    values = [float(r[measure]) for r in items.values()]
                    summaries.append(dict(study=kind, case=case, method=method, metric=measure,
                                          **summarize(values)))
            for base in methods[:-1]:
                values = [float(groups[case,base][b]["seconds"])/float(groups[case,methods[-1]][b]["seconds"])
                          for b in range(1,6)]
                ratios.append(dict(study=kind, case=case, numerator=base, denominator=methods[-1],
                                   **summarize(values)))
    if set(scores)!=set(manifest) or any(len(v)!=1 for v in scores.values()):
        raise ValueError("alignment scores disagree")
    write(opt.out / "summary.csv", sorted(summaries, key=lambda r:(r["study"],r["case"],r["method"],r["metric"])))
    write(opt.out / "paired.csv", sorted(ratios, key=lambda r:(r["study"],r["case"],r["numerator"])))
    for row in data["counters"]:
        if float(row["running_pct"]) < 99 or int(row["value"]) <= 0:
            raise ValueError("invalid hardware counter")


if __name__ == "__main__":
    main()
