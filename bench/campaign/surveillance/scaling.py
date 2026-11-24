#!/usr/bin/env python3
"""Summarize wall time and paired strong scaling of the published matrix executable."""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))
from bench_stats import sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trials", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = pd.read_csv(args.trials)
    measured = data[data.block > 0]
    rows = []
    for algorithm, table in measured.groupby("algorithm"):
        times = table.pivot(index="block", columns="threads", values="wall_s")
        if list(times.columns) != [1, 32] or list(times.index) != list(range(1, 6)) or times.isna().any().any():
            raise ValueError(f"incomplete one-to-32-core comparison for {algorithm}")
        for metric, values in [("wall_1_s", times[1]), ("wall_32_s", times[32]),
                               ("speedup", times[1] / times[32]), ("efficiency", times[1] / times[32] / 32)]:
            stats = sample(values)
            rows.append({"algorithm": algorithm, "metric": metric,
                         **{k: stats[k] for k in ("n", "mean", "sd", "low95", "high95")}})
    pd.DataFrame(rows).to_csv(args.output, index=False, float_format="%.10g")


if __name__ == "__main__":
    main()
