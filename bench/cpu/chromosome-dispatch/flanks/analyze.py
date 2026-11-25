"""Describe paired timings on each chromosome within its warmed process."""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def summarize(trials):
    expected = {("edlib",1), ("wfa2",1), *((v,t) for v in ("reuse","sparse","split") for t in (1,8))}
    cases = {"chr1_MATERNAL":54,"chr1_PATERNAL":111,"chr21_MATERNAL":6,"chr21_PATERNAL":15}
    if set(trials.chrom) != set(cases) or not np.isfinite(trials.seconds).all() or (trials.seconds <= 0).any():
        raise ValueError("Missing case or invalid timing")
    rows = []
    for chrom, records in trials.groupby("chrom", sort=False):
        if set(records.score) != {cases[chrom]}:
            raise ValueError(f"Distance mismatch for {chrom}")
        wide = records.pivot(index="block",columns=["method","threads"],values="seconds")
        if set(wide.index) != set(range(8)) or set(wide.columns) != expected or wide.isna().any().any():
            raise ValueError(f"Incomplete crossover for {chrom}")
        measured = wide.drop(0)
        for (method, threads), times in measured.items():
            reference = measured["reuse", threads]
            ratio = times / reference
            rows.append({"chrom":chrom,"score":cases[chrom],"method":method,"threads":threads,"blocks":len(times),
                         "mean_s":times.mean(),"sd_s":times.std(),"reference":"reuse","reference_threads":threads,
                         "mean_ratio":ratio.mean(),"sd_ratio":ratio.std(),"min_ratio":ratio.min(),"max_ratio":ratio.max()})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trials", type=Path, help="Slurm output directory or combined trials CSV")
    args = parser.parse_args()
    if args.trials.is_dir():
        paths = sorted(args.trials.glob("chr*.csv"))
        if not paths:
            raise ValueError("No trial files found")
        trials = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    else:
        trials = pd.read_csv(args.trials)
    summarize(trials).to_csv(sys.stdout, index=False)
