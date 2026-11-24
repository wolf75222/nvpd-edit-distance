"""Summarize five paired randomized blocks on one and two full EPYC nodes."""
import argparse
from collections import defaultdict
import csv
from pathlib import Path
import statistics


def interval(values):
    values = list(values)
    if len(values) != 5:
        raise ValueError("expected five independent randomized blocks")
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    half = 2.7764451051977987 * sd / 5**0.5
    return mean, sd, mean-half, mean+half


def summarize(directory):
    grouped = defaultdict(dict)
    with (directory / "trials.csv").open() as source:
        for row in csv.DictReader(source):
            key = row["case"], int(row["nodes"]), int(row["height"])
            block = int(row["block"])
            if block in grouped[key]:
                raise ValueError(f"duplicate observation {key}, block {block}")
            grouped[key][block] = row
    if len(grouped) != 15 or any(set(blocks) != set(range(1, 6)) for blocks in grouped.values()):
        raise ValueError("expected 15 configurations and five complete blocks")
    with (directory / "summary.csv").open("w", newline="") as target:
        summary = csv.writer(target)
        summary.writerow(("case","nodes","height","metric","n","mean","sd","low95","high95"))
        for (case,nodes,height), blocks in sorted(grouped.items()):
            for field in ("seconds","gcups","max_rank_rss_kib"):
                summary.writerow((case,nodes,height,field,5,*interval(float(r[field]) for r in blocks.values())))

    with (directory / "scaling.csv").open("w", newline="") as target:
        ratios = csv.writer(target)
        ratios.writerow(("comparison","case","nodes","height","metric","n","mean","sd","low95","high95"))

        def compare(label, case, nodes, height, ref, other, metrics):
            for metric in metrics:
                paired = []
                for block in range(1,6):
                    one = grouped[ref][block]
                    two = grouped[other][block]
                    speed = float(one["seconds"]) / float(two["seconds"])
                    if metric == "strong_efficiency":
                        speed /= 2
                    elif metric == "weak_work_normalized_efficiency":
                        speed *= int(two["cells"]) / (2*int(one["cells"]))
                    paired.append(speed)
                ratios.writerow((label,case,nodes,height,metric,5,*interval(paired)))

        for (case,nodes,height) in sorted(grouped):
            if nodes == 2:
                compare("one_to_two",case,2,height,(case,1,height),(case,2,height),
                        ("strong_speedup","strong_efficiency"))
            if height:
                compare("same_nodes_vs_collective",case,nodes,height,(case,nodes,0),(case,nodes,height),
                        ("speedup_over_collective",))
        compare("384_to_768","matched",2,0,("matched-384",1,0),("matched-768",2,0),
                ("weak_time_ratio","weak_work_normalized_efficiency"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    summarize(parser.parse_args().directory)
