"""Rebuild t(4) summaries of the GH200 calibration and per-allocation paired trials."""

import argparse
import csv
import itertools
import math
from pathlib import Path
from statistics import mean, stdev

from collections import defaultdict

T4 = 2.7764451051977987


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def summarize(values):
    if len(values) != 5 or any(not math.isfinite(v) for v in values):
        raise ValueError("expected five finite process measurements")
    center, sd = mean(values), stdev(values)
    half = T4 * sd / math.sqrt(5)
    return dict(n=5, mean=center, sd=sd, low95=center - half, high95=center + half)


def save(path, records):
    with path.open("w") as stream:
        w = csv.DictWriter(stream, fieldnames=records[0], lineterminator="\n")
        w.writeheader()
        w.writerows(records)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("calibration_job")
    p.add_argument("fusion_job")
    p.add_argument("local_job")
    p.add_argument("classes_job", nargs="?")
    a = p.parse_args()
    rows = read(a.root / a.calibration_job / "calibration.csv")
    curves = defaultdict(list)
    for r in rows:
        k = (r["method"], int(r["steps"]))
        v = float(r["seconds"])
        count = int(r["requested_bytes"] if k[0] == "stream-copy" else r["updates"])
        curves[k].append((r["block"], count / (v * 1e9)))
    if len(rows) != 35 or len(curves) != 7:
        raise ValueError("incomplete calibration")
    cal = []
    for (mode, steps), values in curves.items():
        if {b for b, _ in values} != {str(i) for i in range(1, 6)}:
            raise ValueError("duplicate calibration block")
        cal.append(
            dict(
                job=a.calibration_job,
                method=mode,
                steps=steps,
                metric="GB/s" if mode == "stream-copy" else "Gupdate/s",
                **summarize([v for _, v in values]),
            )
        )
    save(a.root / "calibration-summary.csv", cal)
    summaries = []
    ratios = []
    hashes = {}
    jobs = [
        (a.fusion_job, ("native", "fused", "shuffle", "edlib", "wfa2")),
        (a.local_job, ("native", "fused", "local", "edlib", "wfa2")),
    ]
    if a.classes_job:
        jobs.append((a.classes_job, ("native", "fused", "local", "classes", "edlib")))
    for job, methods in jobs:
        records = read(a.root / job / "trials.csv")
        cfg = defaultdict(dict)
        for r in records:
            case = r["case"]
            block = int(r["block"])
            method = r["method"]
            if r["job"] != job or method not in methods or block in cfg[case, method]:
                raise ValueError("duplicate trial or mismatched job")
            identity = tuple(r[k] for k in ("m", "n", "distance", "input_sha256"))
            if hashes.setdefault(case, identity) != identity:
                raise ValueError("input or score changed")
            cfg[case, method][block] = r
        if len(records) != 200 or len(cfg) != 40:
            raise ValueError("incomplete solver sweep")
        for case in sorted({name for name, _ in cfg}):
            m, n, _, sha = hashes[case]
            for block in range(1, 6):
                if {int(cfg[case, method][block]["order"]) for method in methods} != set(range(1, 6)):
                    raise ValueError("invalid paired block")
            for method in methods:
                obs = [cfg[case, method][i] for i in range(1, 6)]
                times = [float(r["seconds"]) for r in obs]
                rss = [int(r["peak_rss_kib"]) for r in obs]
                metrics = {"seconds": times, "device_bytes": [int(r["device_bytes"]) for r in obs]}
                if method not in ("edlib", "wfa2"):
                    metrics["logical_GCUPS"] = [int(m) * int(n) / v / 1e9 for v in times]
                if min(rss) > 0:
                    metrics["peak_rss_kib"] = rss
                for metric, values in metrics.items():
                    summaries.append(
                        dict(
                            job=job,
                            case=case,
                            input_sha256=sha,
                            method=method,
                            metric=metric,
                            **summarize(values),
                        )
                    )
            gpu_methods = [m for m in methods if m not in ("edlib", "wfa2")]
            for base, test in itertools.combinations(gpu_methods, 2):
                ratios.append(
                    dict(
                        job=job,
                        case=case,
                        input_sha256=sha,
                        numerator=base,
                        denominator=test,
                        **summarize(
                            [
                                float(cfg[case, base][i]["seconds"])
                                / float(cfg[case, test][i]["seconds"])
                                for i in range(1, 6)
                            ]
                        ),
                    )
                )
    save(a.root / "summary.csv", summaries)
    save(a.root / "paired.csv", ratios)


if __name__ == "__main__":
    main()
