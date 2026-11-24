"""Summarize five paired matrix timings after the warm-up block."""
import csv
import math
from pathlib import Path
import statistics
import sys


def estimate(values):
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    return mean, sd, 2.7764451052 * sd / math.sqrt(5)


def main():
    trials = Path(sys.argv[1])
    with trials.open() as stream:
        rows = list(csv.DictReader(stream))
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(("mode", "method", "blocks", "mean_s", "sd_s", "ci95_half_s",
                     "biwfa_over_method", "ci95_half_biwfa_over_method",
                     "split_over_method", "ci95_half_split_over_method"))
    for mode in ("score", "trace"):
        samples = {}
        for method in ("biwfa", "trim-select", "split", "split-filter"):
            group = [r for r in rows if (r["mode"], r["method"]) == (mode, method)]
            if sorted(int(r["block"]) for r in group) != list(range(6)):
                raise ValueError(f"incomplete or repeated blocks for {method}/{mode}")
            samples[method] = [float(r["seconds"]) for r in sorted(group, key=lambda r: int(r["block"]))[1:]]
        for method, times in samples.items():
            ratios = []
            for reference in ("biwfa", "split"):
                mean, _, half = estimate([a / b for a, b in zip(samples[reference], times)])
                ratios.extend((mean, half))
            writer.writerow((mode, method, len(times), *estimate(times), *ratios))


if __name__ == "__main__":
    main()
