#!/usr/bin/env python3
"""Summarize group distances and supported splits from a labeled matrix and an optional public tree."""
import argparse
import csv
import itertools
import math
import statistics
from collections import defaultdict
from pathlib import Path

from phylo_fit import path_matrix, write_table
from restrict_newick import retained_tree
from tree import build_tree, read_matrix


def group_rows(labels, groups, edges, unrooted):
    all_tips = frozenset(range(len(labels)))
    covers = {all_tips}  # Fallback when the requested members span the root.
    for tips, length in edges:
        if length > 1e-9:
            covers.add(tips)
            if unrooted:
                covers.add(all_tips - tips)
    for name in sorted(set(groups.values())):
        members = {i for i, label in enumerate(labels) if groups[label] == name}
        if len(members) < 2:
            yield name, len(members), "", "", ""
        else:
            size = min(len(c) for c in covers if members <= c)
            yield name, len(members), int(size == len(members) and size < len(labels)), size, len(members)/size


def distance_stats(values):
    return len(values), min(values), statistics.median(values), statistics.mean(values), max(values)


def split_set(edges, n):
    universe = frozenset(range(n))
    return {frozenset((tips, universe - tips)) for tips, weight in edges
            if weight > 1e-9 and 1 < len(tips) < n-1}


def compare_tree(public_tree, public_out, labels, groups, d):
    aliases = {s.encode(): s for s in labels}
    root, matched, total = retained_tree(public_tree, aliases)
    if total != len(matched):
        raise ValueError("comparison tree must already be restricted to the manifest")
    shared = [s for s in labels if s in matched]
    if len(shared) < 2:
        raise ValueError("comparison requires at least two shared leaves")
    index = [labels.index(s) for s in shared]
    subset = [[d[i][j] for j in index] for i in index]
    path, edges = path_matrix(root, shared)
    public_out.mkdir(parents=True, exist_ok=True)
    with (public_out / "patristic.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["id", *shared])
        for s, row in zip(shared, path):
            w.writerow([s, *(format(x, ".12g") for x in row)])
    write_table(public_out / "coverage.csv", "id,group,matched",
                ((s, groups[s], int(s in matched)) for s in labels))
    pairs = list(itertools.combinations(range(len(shared)), 2))
    x = [subset[i][j] for i, j in pairs]
    y = [path[i][j] for i, j in pairs]
    def centered(values):
        avg = statistics.mean(values)
        return [v-avg for v in values]
    a, b = centered(x), centered(y)
    denominator = math.sqrt(sum(v*v for v in a)*sum(w*w for w in b))
    coefficient = sum(v*w for v, w in zip(a,b)) / denominator if denominator else ""
    write_table(public_out / "overlap.csv", "shared_tips,pairs,mean_edit_distance,mean_public_path,pearson",
                [(len(shared), len(pairs), statistics.mean(x), statistics.mean(y), coefficient)])
    observed = split_set(edges, len(shared))
    split_rows = []
    group_reports = [("public", *row) for row in group_rows(shared, groups, edges, False)]
    for method in ("upgma", "nj"):
        inferred = build_tree(shared, subset, method)
        _, fit_edges = path_matrix(inferred, shared)
        fitted = split_set(fit_edges, len(shared))
        split_rows.append((method, len(shared), len(fitted), len(observed), len(fitted & observed)))
        group_reports.extend((method, *row) for row in group_rows(shared, groups, fit_edges, method == "nj"))
    write_table(public_out / "shared_splits.csv", "tree,tips,edit_tree_splits,public_splits,shared_splits", split_rows)
    write_table(public_out / "group_clades.csv", "tree,group,members,positive_edge_match,cover_size,purity",
                group_reports)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("matrix", type=Path)
    ap.add_argument("manifest", type=Path, help="TSV with id and group")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--reference-group", help="exclude this label from within/between pooled statistics")
    ap.add_argument("--comparison-tree", type=Path, help="pruned integer-branch Newick, optionally gzipped")
    ap.add_argument("--comparison-out", type=Path)
    args = ap.parse_args()
    if bool(args.comparison_tree) != bool(args.comparison_out):
        ap.error("--comparison-tree and --comparison-out must be specified together")
    labels, d = read_matrix(args.matrix)
    with args.manifest.open() as f:
        records = list(csv.DictReader(f, delimiter="\t"))
    groups = {r["id"]: r["group"] for r in records}
    if len(groups) != len(records) or set(labels) != groups.keys() or not all(groups.values()):
        raise ValueError("each matrix label needs one assigned group")
    args.out.mkdir(parents=True, exist_ok=True)
    grouped, pooled = defaultdict(list), defaultdict(list)
    for i, j in itertools.combinations(range(len(labels)), 2):
        g, h = sorted((groups[labels[i]], groups[labels[j]]))
        grouped[g, h].append(d[i][j])
        if args.reference_group not in (g, h):
            pooled["within" if g == h else "between"].append(d[i][j])
    write_table(args.out / "group_distances.csv", "group_a,group_b,pairs,min,median,mean,max",
                ((*k, *distance_stats(v)) for k, v in sorted(grouped.items())))
    write_table(args.out / "pooled_distances.csv", "comparison,pairs,min,median,mean,max",
                ((k, *distance_stats(v)) for k, v in sorted(pooled.items())))
    reports = []
    for method in ("upgma", "nj"):
        root = build_tree(labels, d, method)
        (args.out / (method + ".nwk")).write_text(root.newick() + ";\n")
        _, edges = path_matrix(root, labels)
        reports.extend((method, *row) for row in group_rows(labels, groups, edges, method == "nj"))
    write_table(args.out / "group_clades.csv", "tree,group,members,positive_edge_match,cover_size,purity", reports)
    if args.comparison_tree:
        compare_tree(args.comparison_tree, args.comparison_out, labels, groups, d)


if __name__ == "__main__":
    main()
