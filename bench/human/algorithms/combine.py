#!/usr/bin/env python3
"""Join allocation-level summaries while preserving their distinct CPU protocols."""
import csv
from pathlib import Path

ROOT = Path("results/human/algorithms")


def main():
    rows = []
    for path in sorted(ROOT.glob("**/summary.csv")):
        if path.parent == ROOT:
            continue
        protocol = "one_core_per_size" if path.parent.parent.name == "arm" else "single_core"
        with path.open() as stream:
            for entry in csv.DictReader(stream):
                entry.setdefault("incomplete", "0")
                rows.append(dict(job=path.parent.name, protocol=protocol, **entry))
    if not rows:
        raise SystemExit("no allocation summaries")
    with (ROOT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
