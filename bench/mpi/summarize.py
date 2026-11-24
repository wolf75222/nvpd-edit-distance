"""Summarize the two seven-block ROMEO experiments, pairing ratios by block."""
import argparse
import csv
from collections import defaultdict
from math import sqrt
from pathlib import Path
from statistics import mean, stdev

# Two-sided 95% interval for seven independent launch blocks (six degrees of freedom).
T_CRITICAL = 2.446911851144969


def summarize(directory):
    groups = defaultdict(dict)
    paired = defaultdict(dict)
    with (directory / "trials.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            case = row["case"]
            nodes, block = int(row["nodes"]), int(row["block"])
            seconds = float(row["seconds"])
            if directory.name == "collectives":
                study = "instrumented" if int(row["instrumented"]) else "control"
                if study == "instrumented":
                    collective = float(row["collective_seconds"])
                    groups[study, case, nodes, "collective_fraction"][block] = collective / seconds
                    groups[study, case, nodes, "collective_seconds"][block] = collective
                paired[case, nodes, block][study] = seconds
            else:
                study = row["study"]
                paired[study, case, block][nodes] = seconds
                groups[study, case, nodes, "equivalent_gcups"][block] = int(row["cells"]) / seconds / 1e9
            leaf = groups[study, case, nodes, "seconds"]
            if block in leaf or seconds <= 0:
                raise ValueError(f"Repeated block or invalid time for {case}")
            leaf[block] = seconds

    for (first, second, block), trial in paired.items():
        if directory.name == "collectives":
            groups["instrumented_vs_control", first, second, "time_ratio"][block] = trial["instrumented"] / trial["control"]
        else:
            ratio = trial[1] / trial[2]
            metric = "speedup" if first == "strong" else "weak_efficiency"
            groups[first, second, 2, metric][block] = ratio
            if first == "strong":
                groups[first, second, 2, "parallel_efficiency"][block] = ratio / 2

    with (directory / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("study", "case", "nodes", "metric", "n", "mean", "sd", "low95", "high95"))
        for key, blocks in sorted(groups.items()):
            if set(blocks) != set(range(1, 8)):
                raise ValueError(f"Incomplete launch blocks for {key}")
            values = list(blocks.values())
            avg, sd = mean(values), stdev(values)
            half = T_CRITICAL * sd / sqrt(len(values))
            writer.writerow((*key, len(values), *(f"{x:.12g}" for x in (avg, sd, avg - half, avg + half))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing trials.csv")
    summarize(parser.parse_args().directory)
