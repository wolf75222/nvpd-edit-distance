#!/usr/bin/env python3
"""Measure three-/four-point gaps and fit two distance trees to a square matrix."""
import argparse
import csv
import itertools
import math
import statistics
from collections import Counter
from pathlib import Path
from tree import build_tree, read_matrix


def write_table(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields.split(","))
        for row in rows:
            writer.writerow(format(v, ".12g") if isinstance(v, float) else v for v in row)


def path_matrix(root, labels):
    position = {name: i for i, name in enumerate(labels)}
    n = len(labels)
    path = [[0.0] * n for _ in labels]
    edges = []

    def visit(node):
        if not node.children:
            return {position[node.label]: 0.0}
        seen = {}
        for child, length in node.children:
            new = {i: d + length for i, d in visit(child).items()}
            for i, old_d in seen.items():
                for j, new_d in new.items():
                    path[i][j] = path[j][i] = old_d + new_d
            edges.append((frozenset(new), length))
            seen.update(new)
        return seen
    if len(visit(root)) != n:
        raise ValueError("tree leaves differ from matrix")
    return path, edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    labels, d = read_matrix(args.matrix)
    if len(labels) < 2:
        parser.error("tree fit requires at least two leaves")
    rows = range(len(labels))
    pairs = list(itertools.combinations(rows, 2))
    # A zero-distance class has identical distance rows in a pseudometric.
    reps = []
    for i in rows:
        duplicates = [j for j in reps if d[i][j] == 0]
        if duplicates and d[i] != d[duplicates[0]]:
            raise ValueError("zero-distance leaves have different distance rows")
        if not duplicates:
            reps.append(i)
    keep = set(reps)
    hist = {"ultrametric": Counter(), "four_point": Counter()}
    collapsed = {"ultrametric": Counter(), "four_point": Counter()}
    for a, b, c in itertools.combinations(rows, 3):
        _, x, y = sorted((d[a][b], d[a][c], d[b][c]))
        hist["ultrametric"][y - x] += 1
        if a in keep and b in keep and c in keep:
            collapsed["ultrametric"][y - x] += 1
    for a, b, c, e in itertools.combinations(rows, 4):
        _, x, y = sorted((d[a][b] + d[c][e], d[a][c] + d[b][e], d[a][e] + d[b][c]))
        hist["four_point"][y - x] += 1
        if a in keep and b in keep and c in keep and e in keep:
            collapsed["four_point"][y - x] += 1
    write_table(args.out / "gaps.csv", "criterion,gap,count",
                ((key, gap, count) for key, counts in hist.items() for gap, count in sorted(counts.items())))
    write_table(args.out / "duplicate_sensitivity.csv", "sample,tips,criterion,combinations,violations,max_gap",
                ((sample, n, name, sum(counts.values()), sum(counts.values())-counts[0], max(counts))
                 for sample, n, tests in (("all", len(labels), hist), ("distinct", len(keep), collapsed))
                 for name, counts in tests.items() if counts))
    fits, branches, pred = [], [], {}
    for method in ("upgma", "nj"):
        path, edges = path_matrix(build_tree(labels, d, method), labels)
        pred[method] = path
        residuals = [path[i][j] - d[i][j] for i, j in pairs]
        fits.append((method, len(pairs), sum(abs(r) < 1e-9 for r in residuals),
                     statistics.mean(residuals), statistics.mean(map(abs, residuals)),
                     math.sqrt(statistics.mean(r*r for r in residuals)), max(map(abs, residuals))))
        for e, (tips, length) in enumerate(edges):
            vals = [r for (i, j), r in zip(pairs, residuals) if (i in tips) != (j in tips)]
            branches.append((method, e, len(tips), length, len(vals), statistics.mean(vals),
                             math.sqrt(statistics.mean(r*r for r in vals)),
                             "|".join(labels[i] for i in sorted(tips))))
    write_table(args.out / "fit.csv", "tree,pairs,zero_residual,mean_residual,mae,rmse,max_abs", fits)
    write_table(args.out / "branches.csv", "tree,edge,descendants,length,crossing_pairs,mean_residual,rmse,tip_ids",
                branches)
    write_table(args.out / "pair_fit.csv", "id_a,id_b,edit_distance,upgma_path,nj_path",
                ((labels[i], labels[j], d[i][j], pred["upgma"][i][j], pred["nj"][i][j]) for i,j in pairs))


if __name__ == "__main__":
    main()
