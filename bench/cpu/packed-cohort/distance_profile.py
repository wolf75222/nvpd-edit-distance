"""Measure the exact edit-score distribution of the nested weak-scaling workloads."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    scores = pd.read_csv(args.matrix, sep="\t", index_col=0)
    if not scores.index.equals(scores.columns) or not scores.equals(scores.T):
        raise ValueError("asymmetric or mislabeled matrix")
    values = scores.to_numpy()
    if np.any(values < 0) or np.any(values.diagonal()):
        raise ValueError("invalid distance matrix")
    observed = np.bincount(values[np.triu_indices(len(scores), 1)])
    expected = pd.read_csv(args.results / "distances.csv").set_index("distance").pairs
    if not np.array_equal(observed, expected.reindex(range(len(observed)), fill_value=0)):
        raise ValueError("matrix disagrees with the ROMEO Edlib histogram")
    designs = sorted(args.results.glob("weak-[0-9]*.csv")) + [args.results / "trials-64.csv"]
    rows = []
    for design in designs:
        case = pd.read_csv(design).iloc[0]
        n = int(case.input_rows)
        d = values[np.triu_indices(n, 1)]
        rows.append((int(case.threads), n, d.size, d.mean(), (d*d).mean(), np.median(d)))
    stats = pd.DataFrame(rows, columns=["threads", "haplotypes", "pairs", "mean_distance", "mean_squared_distance", "median_distance"])
    stats.sort_values("threads").to_csv(args.results/"weak-distance-profile.csv", index=False, float_format="%.8f")


if __name__ == "__main__":
    main()
