"""Run eight randomized process blocks; each process warms once and measures ten calls."""
import argparse
import csv
import hashlib
import math
from pathlib import Path
import random
import statistics
import subprocess

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("bin", "inputs", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--astar", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    programs = {k: args.bin / v for k, v in (
        ("quicked-stock", "measure-quicked-stock"), ("quicked-certified", "measure-quicked"),
        ("edlib", "measure-quicked-native"), ("wfa2", "measure-quicked-native"))}
    if args.astar:
        programs.update({f"astarpa2-{x}": args.astar for x in ("simple", "full")})
    manifest = (args.inputs / "cases.csv").read_bytes()
    (args.out / "cases.csv").write_bytes(manifest)
    cases = list(csv.DictReader(manifest.decode().splitlines()))
    for case in cases:
        path = args.inputs / (case["case"] + ".fa")
        if hashlib.sha256(path.read_bytes()).hexdigest() != case["sha256"]:
            raise ValueError(path)
        oracle = subprocess.check_output([programs["edlib"], path, "edlib", "1"], text=True)
        if int(oracle.split(",")[1]) != int(case["exact"]):
            raise ValueError(f"Independent distance changed: {case}")
    cells = [(c, v) for c in cases for v in programs]
    order = random.Random(216814)
    groups = {}
    with (args.out / "trials.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow("case variant block call seconds score exact correct source retries".split())
        for block in range(1, 9):
            order.shuffle(cells)
            for c, variant in cells:
                name, exact = c["case"], int(c["exact"])
                stdout = subprocess.check_output([programs[variant], args.inputs / (name + ".fa"),
                                                  variant, "10"], text=True, timeout=120)
                observations = list(csv.reader(stdout.splitlines()))
                if len(observations) != 10:
                    raise ValueError((name, variant, stdout))
                times, scores = [], []
                for expected, (call, score, seconds, source, retries) in enumerate(observations, 1):
                    ok = int(score) == exact
                    if int(call) != expected or (variant != "quicked-stock" and not ok) or float(seconds) <= 0:
                        raise ValueError((name, variant, call, score, source))
                    writer.writerow((name, variant, block, call, seconds, score, exact, int(ok), source, retries))
                    times.append(float(seconds))
                    scores.append(int(score))
                group = groups.setdefault((name, variant), ([], [], exact))
                group[0].append(statistics.mean(times))
                group[1].extend(scores)
                f.flush()
    with (args.out / "summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow("case variant blocks calls mean_seconds sd_seconds low95_seconds high95_seconds exact observed correct_calls".split())
        for (name, variant), (means, scores, exact) in sorted(groups.items()):
            mean, sd = statistics.mean(means), statistics.stdev(means)
            half = 2.364624251 * sd / math.sqrt(len(means))  # t(.975,7)
            writer.writerow([name, variant, len(means), len(scores), mean, sd, mean-half, mean+half,
                             exact, " ".join(map(str, sorted(set(scores)))), scores.count(exact)])

if __name__ == "__main__":
    main()
