#!/usr/bin/env python3
"""Build a UPGMA or neighbor-joining Newick tree from a labeled distance TSV."""

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


def read_matrix(path):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader, [])
        labels = header[1:]
        if header[:1] != ["id"] or not labels or "" in labels or len(set(labels)) != len(labels):
            raise ValueError("the header must contain id followed by unique, nonempty labels")
        rows = {}
        for row in reader:
            if len(row) != len(labels) + 1 or row[0] in rows:
                raise ValueError("wrong row width or duplicate row label")
            rows[row[0]] = list(map(float, row[1:]))
    if set(rows) != set(labels):
        raise ValueError("row and column labels differ")
    matrix = [rows[label] for label in labels]
    for i, row in enumerate(matrix):
        if row[i] != 0 or any(not math.isfinite(d) or d < 0 for d in row):
            raise ValueError("distances must be finite and nonnegative, with a zero diagonal")
        if any(d != matrix[j][i] for j, d in enumerate(row)):
            raise ValueError("matrix must be symmetric")
    return labels, matrix


@dataclass
class Node:
    label: str = ""
    children: tuple = ()
    size: int = 1
    height: float = 0.0

    def newick(self):
        if self.children:
            branches = (f"{node.newick()}:{length:.12g}" for node, length in self.children)
            return "(" + ",".join(branches) + ")"
        if re.fullmatch(r"[a-zA-Z0-9.-]+", self.label):
            return self.label
        return "'" + self.label.replace("'", "''") + "'"


def build_tree(labels, matrix, method):
    active = {i: Node(label) for i, label in enumerate(labels)}
    distances = [dict(enumerate(row)) for row in matrix]
    while len(active) > 1:
        pairs = combinations(active, 2)
        if method == "nj" and len(active) > 2:
            divisor = len(active) - 2
            totals = {i: sum(distances[i].values()) for i in active}
            i, j = min(pairs, key=lambda p: divisor * distances[p[0]][p[1]] - totals[p[0]] - totals[p[1]])
            limb_i = (distances[i][j] + (totals[i] - totals[j]) / divisor) / 2
        else:
            i, j = min(pairs, key=lambda p: distances[p[0]][p[1]])
            limb_i = distances[i][j] / 2
        dij = distances[i][j]
        left, right = active.pop(i), active.pop(j)
        size = left.size + right.size
        if method == "upgma":
            height = dij / 2
            children = ((left, height - left.height), (right, height - right.height))
        else:
            height = 0
            children = ((left, limb_i), (right, dij - limb_i))
        parent = len(distances)
        new_row = {parent: 0.0}
        for k in active:
            di, dj = distances[k].pop(i), distances[k].pop(j)
            if method == "upgma":
                d = (left.size * di + right.size * dj) / size
            else:
                d = (di + dj - dij) / 2
            distances[k][parent] = new_row[k] = d
        distances[i].clear()
        distances[j].clear()
        distances.append(new_row)
        active[parent] = Node(children=children, size=size, height=height)
    return next(iter(active.values()))


def report_groups(path, labels, root, method):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)
        if not {"id", "group"} <= set(reader.fieldnames or []):
            raise ValueError("group metadata needs columns id and group")
        groups, seen = {}, set()
        for row in reader:
            name, group = row["id"], row["group"]
            if not name or not group or name in seen:
                raise ValueError("expected one nonempty group per unique ID")
            seen.add(name)
            groups.setdefault(group, set()).add(name)
    if seen != set(labels):
        raise ValueError("group IDs must match the matrix labels")

    clades = set()

    def descend(node):
        tips = set() if node.children else {node.label}
        for child, length in node.children:
            members = descend(child)
            if length > 0:
                clades.add(frozenset(members))
                if method == "nj":
                    clades.add(frozenset(seen - members))
            tips.update(members)
        return tips

    descend(root)
    print("group\tmembers\tpositive_edge_match", file=sys.stderr)
    for group, members in sorted(groups.items()):
        print(group, len(members), int(members in clades), sep="\t", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--method", choices=("upgma", "nj"), default="upgma")
    parser.add_argument("--out", type=Path, help="Newick output (default: stdout)")
    parser.add_argument("--groups", type=Path, help="TSV with id and group; report matching edges on stderr")
    args = parser.parse_args()
    try:
        labels, matrix = read_matrix(args.matrix)
        root = build_tree(labels, matrix, args.method)
        if args.groups:
            report_groups(args.groups, labels, root, args.method)
        newick = root.newick() + ";\n"
        if args.out:
            args.out.write_text(newick, encoding="utf-8")
        else:
            sys.stdout.write(newick)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"{parser.prog}: {exc}\n")


if __name__ == "__main__":
    main()
