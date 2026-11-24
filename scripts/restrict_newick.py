#!/usr/bin/env python3
"""Prune a plain or gzipped accession-labeled Newick tree with integer branch lengths."""
import argparse
import csv
import gzip
import re
from pathlib import Path
from tree import Node


def retained_tree(path, aliases):
    # Keep only ancestors of selected leaves. Other branches need no node objects.
    frames, found = [], set()
    pending, closed, root = None, False, None
    total = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as stream:
        tail = b""
        for chunk in iter(lambda: stream.read(2**20), b""):
            data = tail + chunk
            begin = 0
            for sep in re.finditer(rb"[(),;]", data):
                token = data[begin:sep.start()].strip()
                begin = sep.end()
                delim = sep[0]
                if delim == b"(":
                    if token or closed:
                        raise ValueError("unsupported Newick node prefix")
                    frames.append([])
                    continue
                if not token and closed and delim == b";" and not frames:
                    token = b":0"
                title, colon, scalar = token.rpartition(b":")
                if not colon:
                    raise ValueError("expected explicit branch lengths")
                if not closed:
                    total += 1
                    if title.startswith(b"'") and title.endswith(b"'"):
                        title = title[1:-1].replace(b"''", b"'")
                    hits = {aliases[field] for field in title.split(b"|") if field in aliases}
                    if len(hits) > 1 or hits & found:
                        raise ValueError("ambiguous accession at leaf")
                    found.update(hits)
                    pending = (Node(label=hits.pop()), 0) if hits else None
                if pending is not None:
                    node, length = pending
                    item = (node, length + int(scalar))
                    if frames:
                        frames[-1].append(item)
                    else:
                        root = node
                pending, closed = None, False
                if delim == b")":
                    members = frames.pop()
                    if len(members) == 1:
                        pending = members[0]
                    elif members:
                        pending = (Node(children=tuple(members), size=sum(k.size for k, _ in members)), 0)
                    closed = True
            tail = data[begin:]
        if tail.strip() or frames or root is None:
            raise ValueError("incomplete or empty restricted tree")
    return root, found, total


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("tree", type=Path, help="plain or gzip-compressed uncondensed Newick")
    p.add_argument("manifest", type=Path, help="TSV with id column")
    p.add_argument("--alias", action="append", default=[], help="equivalent source ID=manifest ID")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    with args.manifest.open() as stream:
        ids = [r["id"] for r in csv.DictReader(stream, delimiter="\t")]
    aliases = {v.encode(): v for v in ids}
    for alias in args.alias:
        a, b = alias.split("=")
        if b not in ids:
            raise ValueError("alias target absent from manifest")
        aliases[a.encode()] = b
    tree, found, tips = retained_tree(args.tree, aliases)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(tree.newick() + ";\n")
    print(f"source leaves={tips} matched={len(found)}/{len(ids)}")
    if len(found) < len(ids):
        print("Unmatched accessions:", ", ".join(x for x in ids if x not in found))


if __name__ == "__main__":
    main()
