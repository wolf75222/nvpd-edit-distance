"""Summarize five complete randomized blocks and their within-block ratios."""
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def check_inputs(folder):
    human = Path(__file__).resolve().parents[3] / "results/human"
    expected_hash = {}
    for file, column in ((human/"chromosome-panel/manifest.tsv", "input_sha256"),
                         (human/"whole/manifest.tsv", "sha256")):
        with file.open() as stream:
            for row in csv.DictReader(stream,delimiter="\t"):
                expected_hash[row["case"]] = row[column]
    for line in (folder/"inputs.sha256").read_text().splitlines():
        digest, name = line.split()
        if digest != expected_hash[Path(name).stem]:
            raise ValueError(f"input checksum disagrees with source manifest: {name}")
    scores = {}
    for file in (human/"chromosome-panel/trials.csv",human/"algorithms/838910/trials.csv"):
        with file.open() as stream:
            for row in csv.DictReader(stream):
                if row["method"] not in ("edlib","wfa2"):
                    continue
                if scores.setdefault(row["case"], row["distance"]) != row["distance"]:
                    raise ValueError(f"independent scores disagree for {row['case']}")
    return scores


def summarize(folder):
    with (folder / "trials.csv").open() as f:
        rows = list(csv.DictReader(f))
    groups = defaultdict(dict)
    scores = check_inputs(folder)
    capped = set()
    if (folder/"censored.csv").is_file():
        with (folder/"censored.csv").open() as f:
            capped = {(r["case"],r["mode"],r["method"],int(r["threads"])) for r in csv.DictReader(f)}
    for row in rows:
        case, mode, alg = (row[k] for k in ("case", "mode", "method"))
        threads, block = int(row["threads"]), int(row["block"])
        obs = groups[case, mode, alg, threads]
        if block in obs or scores[case] != row["distance"]:
            raise ValueError(row)
        obs[block] = row

    with (folder / "summary.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case","mode","method","threads","metric","n","mean","sd","low95","high95"))

        def record(key, metric, values):
            if len(values) != 5:
                raise ValueError((key,metric,len(values)))
            mean, sd = statistics.mean(values), statistics.stdev(values)
            half = 2.7764451052 * sd / math.sqrt(5)
            writer.writerow((*key, metric, 5, mean, sd, mean-half, mean+half))

        for key in capped:
            groups.pop(key, None)
        for key, observations in sorted(groups.items()):
            case, mode, alg, threads = key
            if sorted(observations) != list(range(1,6)):
                raise ValueError(f"incomplete configuration {key}")
            for metric in ("seconds","peak_rss_kib"):
                record(key, metric, [float(row[metric]) for row in observations.values()])
            if alg != "split":
                continue
            for baseline, width in (("split",1),("split-all",8),("trim-only",1),("biwfa",1),("wfa2",1),("edlib",1)):
                reference = groups.get((case,mode,baseline,width))
                if reference:
                    ratios = [float(reference[i]["seconds"])/float(obs["seconds"])
                              for i,obs in observations.items()]
                    record(key, f"speedup_{baseline}", ratios)
                    if baseline == "split":
                        record(key, "efficiency", [r/threads for r in ratios])


if __name__ == "__main__":
    summarize(Path(sys.argv[1]))
