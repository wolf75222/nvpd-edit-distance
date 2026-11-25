#!/usr/bin/env python3
"""Export the hottest application loop from native ONEVIEW experiment directories."""
import argparse
import csv
import json
from pathlib import Path
import sys


def records(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter=";"))


def describe(profile):
    run = profile / "shared/run_0"
    loops = records(run / "expert_loops.csv")
    hot = max((r for r in loops if r["Level"] in ("Single", "Innermost")),
              key=lambda r: float(r["Coverage (% app. time)"]))
    module = hot["Module"]
    key = f"{module}_{int(hot['ID'].split()[-1])}"
    details = profile / "static_data"
    assembly = [line["insn"].split()[0].upper()
                for line in records(details / "asm" / f"{key}.csv")]
    cqa = records(details / "cqa" / f"{key}_cqa.csv")[0]
    advice = json.loads((details / "cqa" / f"{key}_text.json").read_text())
    # The shipped kernels have one path in their hot loops.
    paths = advice["_cqa_text_report"]["paths"]
    if len(paths) != 1:
        raise ValueError(f"Multiple machine-code paths in {key}")
    global_metrics = {r["metric"]: r["value"] for r in records(run / "global_metrics.csv")}
    categories = next(r for r in records(run / "lprof_categorization.csv")
                      if r["Node"] == r["Process"] == r["Thread"] == "_ALL")
    omp_percent = 100 * float(categories["OMP"]) / float(categories["#total"])
    return {
        "directory": profile.name, "module": module, "loop": key, "source": hot["Source Location"],
        "event": categories["Metric"], "seconds_profiled": global_metrics["profiled_time"],
        "omp_percent": f"{omp_percent:.4f}", "loop_percent": hot["Coverage (% app. time)"],
        "vectorization_percent": hot["Vectorization Ratio (%)"],
        "width_efficiency_percent": hot["Vectorization Efficiency (%)"],
        "model_cycles_per_asm_iteration": hot["CQA cycles"],
        "instruction_count": cqa["nb instructions"],
        "indirect_streams": cqa["stride_indirect"],
        "gather_instructions": sum("GATHER" in insn for insn in assembly),
        "permutation_instructions": sum("PERM" in insn for insn in assembly),
        "hints": " | ".join(r["title"] for r in paths[0]["hint"]),
        "potential": " | ".join(r["title"] for r in paths[0]["potential"]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", nargs="+", type=Path)
    args = parser.parse_args()
    rows = [describe(path) for path in args.profiles]
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
