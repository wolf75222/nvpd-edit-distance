"""Check widths, scores and placements before aggregating independent process blocks."""

import argparse
import csv
from collections import defaultdict
from itertools import combinations
from math import isfinite, sqrt
from pathlib import Path
from statistics import mean, stdev

from measure import write

METHODS = ("native", "kokkos-rows", "kokkos-opt")
T95_DF4 = 2.7764451051977987


def stats(values, name):
    if len(values) != 5 or any(not isfinite(v) for v in values):
        raise ValueError("expected five finite independent observations")
    mu, sd = mean(values), stdev(values)
    half = T95_DF4 * sd / sqrt(len(values))
    return {
        f"mean_{name}": mu,
        f"sd_{name}": sd,
        f"ci95_low_{name}": mu - half,
        f"ci95_high_{name}": mu + half,
    }


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("epyc", type=Path)
    p.add_argument("gh200", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    refs = read(a.epyc / "references.csv")
    if refs != read(a.gh200 / "references.csv") or len(refs) != 10:
        raise ValueError("the independent reference tables do not agree")
    source = {r["case"]: r for r in refs}
    observations = {}
    cpus = defaultdict(set)
    for machine, folder in (("epyc", a.epyc), ("gh200", a.gh200)):
        for path in sorted(folder.glob("lane-*.csv")):
            rows = read(path)
            if len({x["cpu_ids"] for x in rows}) != 1:
                raise ValueError("worker placement changed")
            affinity = {int(c) for c in rows[0]["cpu_ids"].split(";")}
            if cpus[machine] & affinity or len(affinity) != int(rows[0]["workers"]):
                raise ValueError("overlapping or invalid worker masks")
            cpus[machine].update(affinity)
            for x in rows:
                case = x["case"]
                alg = x["algorithm"]
                w = int(x["workers"])
                b = int(x["block"])
                key = (machine, case, w, alg, b)
                if (
                    key in observations
                    or x["target"] != machine
                    or x["execution_space"] != ("OpenMP" if machine == "epyc" else "Cuda")
                ):
                    raise ValueError(f"invalid trial identity at {key}")
                if any(x[k] != source[case][k] for k in ("m", "n", "region", "distance", "input_sha256")):
                    raise ValueError("oracle mismatch")
                if not float(x["seconds"]) > 0 or int(x["peak_rss_kib"]) <= 0:
                    raise ValueError("invalid process measurement")
                m, n = int(x["m"]), int(x["n"])
                explicit = 0 if machine == "epyc" else 6 * 4 * (max(m, n) + 1)
                if int(x["device_bytes"]) != explicit or int(x["cells"]) != m * n:
                    raise ValueError("incorrect device allocation or cell count")
                observations[key] = x
    expected = {
        (machine, r["case"], w, method, b)
        for r in refs
        for machine in ("epyc", "gh200")
        if machine == "gh200" or int(r["m"]) <= 65536
        for w in ((1, 4, 8) if machine == "epyc" else (1,))
        for method in METHODS
        for b in range(1, 6)
    }
    if set(observations) != expected:
        raise ValueError("incomplete workload coverage")
    trials = []
    summary = []
    paired = []
    scaling = []
    for machine in ("epyc", "gh200"):
        for ref in refs:
            name = ref["case"]
            basetimes = {}
            if machine == "epyc" and int(ref["m"]) > 65536:
                continue
            for w in (1, 4, 8) if machine == "epyc" else (1,):
                common = dict(
                    target=machine,
                    case=name,
                    region=ref["region"],
                    m=ref["m"],
                    n=ref["n"],
                    distance=ref["distance"],
                    input_sha256=ref["input_sha256"],
                    workers=w,
                    blocks=5,
                )
                times = {}
                for method in METHODS:
                    runs = [observations[machine, name, w, method, b] for b in range(1, 6)]
                    trials.extend(runs)
                    times[method] = [float(x["seconds"]) for x in runs]
                    means = stats(times[method], "seconds")
                    summary.append(
                        dict(
                            **common,
                            algorithm=method,
                            **means,
                            **stats([float(x["gcups"]) for x in runs], "gcups"),
                            **stats([float(x["peak_rss_kib"]) for x in runs], "peak_rss_kib"),
                            device_bytes=runs[0]["device_bytes"],
                            host_workspace_array_bytes=runs[0]["host_workspace_array_bytes"],
                        )
                    )
                    if machine == "epyc":
                        if w == 1:
                            basetimes[method] = means["mean_seconds"]
                        speedup = basetimes[method] / means["mean_seconds"]
                        scaling.append(
                            dict(
                                **common,
                                algorithm=method,
                                speedup_vs_one_worker=speedup,
                                efficiency=speedup / w,
                            )
                        )
                for x, y in combinations(METHODS, 2):
                    paired.append(
                        dict(
                            **common,
                            numerator=x,
                            denominator=y,
                            **stats([u / v for u, v in zip(times[x], times[y])], "ratio"),
                        )
                    )
    a.out.mkdir(parents=True, exist_ok=True)
    for fname, rows in [
        ("references", refs),
        ("trials", trials),
        ("summary", summary),
        ("paired", paired),
        ("scaling", scaling),
    ]:
        write(a.out / (fname + ".csv"), rows)
    print(
        f"{len(trials)} fresh processes, {len(summary)} configurations, {len(paired)} paired comparisons."
    )


if __name__ == "__main__":
    main()
