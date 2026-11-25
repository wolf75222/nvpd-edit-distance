"""Summarize the two seven-block ROMEO experiments, pairing ratios by block."""
import argparse
import csv
from collections import defaultdict
from math import isfinite, sqrt
from pathlib import Path
from statistics import mean, stdev
from tempfile import NamedTemporaryFile

# Two-sided 95% interval for seven independent launch blocks (six degrees of freedom).
T_CRITICAL = 2.446911851144969


def summarize(directory):
    groups = defaultdict(dict)
    paired = defaultdict(dict)
    identities, cells_per_node, jobs = {}, {}, set()
    with (directory / "trials.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or ())
        collectives = "instrumented" in fields
        if collectives == ("study" in fields):
            raise ValueError("Expected collective or scaling trial columns")
        for row in reader:
            case = row["case"]
            nodes, block = int(row["nodes"]), int(row["block"])
            seconds = float(row["seconds"])
            jobs.add(row["job"])
            if (not isfinite(seconds) or seconds <= 0 or nodes not in (1, 2)
                    or int(row["physical_cores"]) != 192 * nodes):
                raise ValueError(f"Invalid time or node allocation for {case}")
            if collectives:
                if row["instrumented"] not in ("0", "1"):
                    raise ValueError(f"Invalid instrumentation flag for {case}")
                study = "instrumented" if int(row["instrumented"]) else "control"
                key = ("collectives", case)
                identity = (row["input_sha256"], row["distance"])
                if study == "instrumented":
                    collective = float(row["collective_seconds"])
                    if not isfinite(collective) or collective < 0:
                        raise ValueError(f"Invalid collective time for {case}")
                    groups[study, case, nodes, "collective_fraction"][block] = collective / seconds
                    groups[study, case, nodes, "collective_seconds"][block] = collective
                paired[case, nodes, block][study] = seconds
            else:
                study = row["study"]
                cells = int(row["cells"])
                if study not in ("strong", "weak") or cells <= 0:
                    raise ValueError(f"Invalid scaling workload for {case}")
                key = (study, case, nodes if study == "weak" else None)
                identity = (row["input_sha256"], row["oracle"], cells, int(row["pairs"]))
                if study == "weak":
                    reference_cells, reference_nodes = cells_per_node.setdefault(case, (cells, nodes))
                    if reference_cells * nodes != cells * reference_nodes:
                        raise ValueError(f"Weak-scaling cells per node changed for {case}")
                paired[study, case, block][nodes] = seconds
                groups[study, case, nodes, "equivalent_gcups"][block] = cells / seconds / 1e9
            if identities.setdefault(key, identity) != identity:
                raise ValueError(f"Input or reference result changed for {case}")
            leaf = groups[study, case, nodes, "seconds"]
            if block in leaf:
                raise ValueError(f"Repeated block for {case}")
            leaf[block] = seconds

    if len(jobs) != 1:
        raise ValueError("Expected trials from one Slurm allocation")
    for (first, second, block), trial in paired.items():
        if trial.keys() != ({"control", "instrumented"} if collectives else {1, 2}):
            raise ValueError(f"Incomplete comparison for {first}/{second}, block {block}")
        if collectives:
            groups["instrumented_vs_control", first, second, "time_ratio"][block] = trial["instrumented"] / trial["control"]
        else:
            ratio = trial[1] / trial[2]
            metric = "speedup" if first == "strong" else "weak_efficiency"
            groups[first, second, 2, metric][block] = ratio
            if first == "strong":
                groups[first, second, 2, "parallel_efficiency"][block] = ratio / 2

    stream = NamedTemporaryFile("w", dir=directory, newline="", delete=False)
    temporary = Path(stream.name)
    try:
        with stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(("study", "case", "nodes", "metric", "n", "mean", "sd", "low95", "high95"))
            for key, blocks in sorted(groups.items()):
                if set(blocks) != set(range(1, 8)):
                    raise ValueError(f"Incomplete launch blocks for {key}")
                values = list(blocks.values())
                avg, sd = mean(values), stdev(values)
                half = T_CRITICAL * sd / sqrt(len(values))
                writer.writerow((*key, len(values), *(f"{x:.12g}" for x in (avg, sd, avg - half, avg + half))))
        temporary.replace(directory / "summary.csv")
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing trials.csv")
    summarize(parser.parse_args().directory)
