"""Aggregate the three matched repeats within each of five independent FASTA pairs."""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "python"))
from bench_stats import sample


def summarize(trials):
    if (len(trials) != 40 or trials.duplicated(["pair", "block", "method"]).any()
            or set(trials.pair) != set(range(1, 6)) or set(trials.block) != set(range(4))
            or set(trials.method) != {"nvpd-omp", "diagonal"}
            or not np.isfinite(trials.seconds).all() or (trials.seconds <= 0).any()):
        raise ValueError("Expected five pairs, four blocks and two positive durations per block")
    if trials.groupby("pair")[["length_a", "length_b", "score"]].nunique().ne(1).any().any():
        raise ValueError("Scores or input lengths differ")
    paired = trials[trials.block > 0].pivot(index=["pair","block"], columns="method", values="seconds")
    paired["speedup"] = paired.diagonal / paired["nvpd-omp"]
    by_pair = paired.groupby("pair").mean()
    return pd.DataFrame([{"metric": metric, **sample(by_pair[metric])} for metric in by_pair])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trials", type=Path, help="Slurm output directory or combined trials CSV")
    args = parser.parse_args()
    if args.trials.is_dir():
        paths = sorted(args.trials.glob("pair-*.csv"))
        if not paths:
            raise ValueError("No trial files found")
        trials = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    else:
        trials = pd.read_csv(args.trials)
    summarize(trials).to_csv(sys.stdout, index=False)
