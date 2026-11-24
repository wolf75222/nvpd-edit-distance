"""Count logical source-level word requests for interior cells of unfused NvPD."""

import argparse
import csv
import hashlib
from pathlib import Path
from bisect import bisect_right
from array import array


def plan_work(a, b, columns):
    positions = {letter: array("I") for letter in set(a)}
    for j, letter in enumerate(b, 1):
        if letter in positions:
            positions[letter].append(j)
    halo, widths = [], []
    full_work = tail_work = 0
    tail = len(a) % 8
    for start in range(0, len(b) + 1, columns):
        end = min(start + columns - 1, len(b))
        edges = [start]
        for _ in range(7):
            left = edges[-1]
            before = max(0, left - 1)
            for hits in positions.values():
                k = bisect_right(hits, left)
                if k:
                    before = min(before, hits[k - 1] - 1)
            edges.append(before)
        halo.append(start - edges[-1])
        widths.append(end - edges[-1] + 1)
        full_work += sum(end - x + 1 for x in edges)
        tail_work += sum(end - x + 1 for x in edges[:tail])
    cutoff = sorted(widths)[(9 * len(widths) - 1) // 10]
    small = [width for width in widths if width <= cutoff]
    large = [width for width in widths if width > cutoff]
    return (
        len(halo),
        max(halo),
        sum(halo) / len(halo),
        full_work * (len(a) // 8) + tail_work,
        len(small),
        len(large),
        8 * max(small),
        8 * max(large, default=0),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("references", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--trials", type=Path, help="recover the actual local tile widths from allocated bytes"
    )
    args = p.parse_args()
    allocations = {}
    if args.trials:
        with args.trials.open() as f:
            for row in csv.DictReader(f):
                key = (row["case"], row["method"])
                size = int(row["device_bytes"])
                if allocations.setdefault(key, size) != size:
                    raise ValueError("workspace changed")
    with args.references.open() as f, args.out.open("w") as out:
        w = csv.writer(out, lineterminator="\n")
        header = (
            "case",
            "input_sha256",
            "m",
            "n",
            "exact_matches",
            "no_prior_match",
            "other",
            "source_request_bytes",
            "logical_cells_per_requested_byte",
        )
        if allocations:
            header += (
                "output_width",
                "blocks",
                "max_halo",
                "mean_halo",
                "computed_cells_including_boundaries",
                "small_blocks",
                "large_blocks",
                "small_shared_bytes",
                "large_shared_bytes",
            )
        w.writerow(header)
        for r in csv.DictReader(f):
            if int(r["m"]) not in (16384, 65536, 131072, 524288):
                continue
            raw = (args.inputs / (r["case"] + ".fa")).read_bytes()
            if hashlib.sha256(raw).hexdigest() != r["input_sha256"]:
                raise ValueError("FASTA digest changed")
            # The input generator writes one sequence per line.
            h1, a, h2, b = raw.splitlines()
            if not h1.startswith(b">") or not h2.startswith(b">"):
                raise ValueError("invalid reference FASTA")
            if len(a) > len(b):
                a, b = b, a
            cells = len(a) * len(b)
            hit = absent = 0
            for letter in set(a):
                char = bytes([letter])
                count = a.count(char)
                hit += count * b.count(char)
                first = b.find(char)
                absent += count * (len(b) if first < 0 else first)
            other = cells - hit - absent
            # Each cell loads one match index and stores one score. The three
            # branches reference one, two or three previous scores, respectively.
            requested = 4 * (3 * hit + 4 * absent + 5 * other)
            result = (
                r["case"],
                r["input_sha256"],
                len(a),
                len(b),
                hit,
                absent,
                other,
                requested,
                cells / requested,
            )
            if allocations:
                case = r["case"]
                extra = allocations[case, "local"] - allocations[case, "native"]
                widths = [w for w in (256, 1024) if ((len(b) + w) // w) * 8 * 4 == extra]
                if len(widths) != 1:
                    raise ValueError("could not identify the measured CUDA tile width")
                result += (widths[0], *plan_work(a, b, widths[0]))
            w.writerow(result)


if __name__ == "__main__":
    main()
