"""Validate chromosome identities and summarize nine matched blocks per page policy."""
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

CASES = {"chr1": 176, "chr18": 5, "chr20": 30}


def observations(row):
    calls = int(row["calls"])
    core_s = float(row["core_s"])
    mem_s = float(row["mem_s"])
    fills = int(row["all_fills"])
    instructions = int(row["instructions"])
    return {
        "score_ms": 1000 * core_s / calls,
        "map_ms": 1000 * float(row["setup_s"]),
        "thp_fraction": int(row["huge_bytes"]) / int(row["mapped_bytes"]),
        "cycles_per_score": int(row["cycles"]) / calls,
        "instructions_per_score": instructions / calls,
        "ipc": instructions / int(row["cycles"]),
        "l1d_miss_per_score": int(row["l1d_miss"]) / calls,
        "walks_per_score": int(row["dtlb_walks"]) / calls,
        "memory_fill_GBs": 64 * fills / mem_s / 1e9,
        "demand_GBs": 64 * int(row["dem_fills"]) / mem_s / 1e9,
        "hw_prefetch_GBs": 64 * int(row["prefetch_fills"]) / mem_s / 1e9,
    }


def summarize(folder):
    root = Path(__file__).resolve().parents[3]
    with (root / "results/human/chromosome-panel/manifest.tsv").open() as stream:
        source = {r["case"]: r for r in csv.DictReader(stream, delimiter="\t")}
    expected_hashes = {key: source[key]["input_sha256"] for key in CASES}
    observed_hashes = {}
    for line in (folder / "inputs.sha256").read_text().splitlines():
        digest, file = line.split()
        observed_hashes[Path(file).stem] = digest
    if observed_hashes != expected_hashes:
        raise ValueError("FASTA inputs differ from the CHM13 panel")

    data = defaultdict(dict)
    with (folder / "trials.csv").open() as stream:
        for r in csv.DictReader(stream):
            case, solver, page = r["case"], r["solver"], r["pages"]
            block = int(r["block"])
            if not (int(r["distance"]) == CASES[case] and r["m"] == source[case]["m"] and r["n"] == source[case]["n"]):
                raise ValueError(("input or score mismatch", r))
            if min(float(r["core_active"]), float(r["mem_active"])) < .99 or int(r["calls"]) != 3:
                raise ValueError(("incomplete PMU batch", r))
            if block in data[case, solver, page]:
                raise ValueError(("duplicate block", r))
            data[case, solver, page][block] = observations(r)

    blocks = range(1, 10)
    expected = {(c, s, p) for c in CASES for s in ("split", "wfa2") for p in ("base", "thp")}
    if set(data) != expected or any(set(v) != set(blocks) for v in data.values()):
        raise ValueError("incomplete experiment")

    with (folder / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case", "solver", "pages", "metric", "n", "mean", "sd", "low95", "high95"))

        def record(key, metric, values):
            mean = statistics.mean(values)
            sd = statistics.stdev(values)
            # Quantile t(.975, 8), two-sided interval for a mean of nine blocks.
            half = 2.3060041352 * sd / math.sqrt(len(values))
            writer.writerow((*key, metric, len(values), mean, sd, mean-half, mean+half))

        for key, samples in sorted(data.items()):
            case, solver, page = key
            for metric in samples[1]:
                values = [samples[b][metric] for b in blocks]
                record(key, metric, values)
                if page == "thp" and metric in ("score_ms", "map_ms", "walks_per_score", "l1d_miss_per_score", "cycles_per_score"):
                    base = data[case, solver, "base"]
                    record((case, solver, "base/thp"), metric, [base[b][metric] / samples[b][metric] for b in blocks])


def summarize_reuse(folder):
    series = defaultdict(dict)
    with (folder / "reuse.csv").open() as stream:
        for row in csv.DictReader(stream):
            page, k, block = row["pages"], int(row["queries"]), int(row["block"])
            if block in series[page, k] or int(row["distance"]) != CASES["chr1"]:
                raise ValueError(("invalid repeated score", row))
            if float(row["total_s"]) <= float(row["setup_s"]):
                raise ValueError(("invalid timer", row))
            series[page, k][block] = float(row["total_s"])
    if set(series) != {(p, k) for p in ("base", "thp") for k in (1, 192, 384)}:
        raise ValueError("missing chromosome reuse configuration")
    blocks = range(1, 6)
    if any(set(g) != set(blocks) for g in series.values()):
        raise ValueError("incomplete chromosome reuse block")
    with (folder / "reuse-summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("queries", "metric", "n", "mean", "sd", "low95", "high95"))
        for k in (1, 192, 384):
            base = series["base", k]
            thp = series["thp", k]
            metrics = {"base_s": base.values(), "thp_s": thp.values(),
                       "base_minus_thp_ms": [1000 * (base[b] - thp[b]) for b in blocks]}
            for key, values in metrics.items():
                average, sd = statistics.mean(values), statistics.stdev(values)
                half = 2.7764451052 * sd / math.sqrt(5)
                writer.writerow((k, key, 5, average, sd, average-half, average+half))


if __name__ == "__main__":
    output = Path(sys.argv[1])
    summarize(output)
    summarize_reuse(output)
