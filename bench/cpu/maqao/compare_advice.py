#!/usr/bin/env python3
"""Compare ONEVIEW's full CQA report on identical binaries at one and four cores."""
import argparse
import csv
from difflib import unified_diff
from hashlib import sha256
import json
from pathlib import Path

from hotloops import describe, records


def normalize(value):
    """Ignore text spacing and JSON key order, retaining lists and linked detail objects."""
    if isinstance(value, str):
        return " ".join(value.split())
    if isinstance(value, list):
        return [normalize(item) for item in value]
    if isinstance(value, dict):
        return {key: normalize(value[key]) for key in sorted(value)}
    return value


def titles(report):
    levels = ("brief", "gain", "potential", "hint")
    return [{level: [entry["title"] for entry in path[level]] for level in levels} for path in report["paths"]]


def read_profile(path, version, kernel, threads):
    summary = describe(path)
    run = path / "shared/run_0"
    settings = {r["metric"]: r["value"] for r in records(run / "local_vars.csv")}
    metrics = {r["metric"]: r["value"] for r in records(run / "global_metrics.csv")}
    if (settings["maqao_version"], int(metrics["nb_threads"]), summary["module"]) != (version, threads, kernel):
        raise ValueError(f"Export does not match version, thread count or binary in {path}")
    report_file = path / "static_data/cqa" / (summary["loop"] + "_text.json")
    report = normalize(json.loads(report_file.read_text())["_cqa_text_report"])
    canonical = json.dumps(report, ensure_ascii=False, indent=2)
    fields = ("loop", "event", "seconds_profiled", "loop_percent", "model_cycles_per_asm_iteration",
              "instruction_count", "gather_instructions", "permutation_instructions")
    binary = (path / "binaries" / kernel).read_bytes()
    row = {"version": version, "build": settings["maqao_build"], "kernel": kernel, "threads": threads,
           "architecture": settings["uarchitecture"], "binary_sha256": sha256(binary).hexdigest(),
           **{key: summary[key] for key in fields}, "cells_per_asm_iteration": 8 if kernel == "original" else 16,
           "cqa_report_sha256": sha256(canonical.encode()).hexdigest()}
    return row, report, canonical


def object_inventory(profile):
    exported = {}
    for path in (profile / "static_data/cqa").glob("*_text.json"):
        key = path.name.removesuffix("_text.json")
        disassembly = records(profile / "static_data/asm" / f"{key}.csv")
        code = [(r["address"], " ".join(r["insn"].split())) for r in disassembly]
        if not code:
            raise ValueError(f"Missing machine instructions for {key}")
        payload = normalize(json.loads(path.read_text())["_cqa_text_report"])
        exported[key] = (sha256(json.dumps(code).encode()).hexdigest(),
                         json.dumps(payload, ensure_ascii=False, indent=2))
    sampled = {}
    for kind in ("loop", "function"):
        files = list((profile / "shared").glob(f"lprof*run_0/lprof_{kind}s.csv"))
        if len(files) != 1:
            raise ValueError(f"Expected one aggregate LProf {kind} CSV in {profile}")
        for row in records(files[0]):
            prefix = "fct_" if kind == "function" else ""
            key = f"{prefix}{row['module']}_{row[kind + '_id']}"
            sampled[key] = row["inc_time_p"]
    return exported, sampled


def compare(root, versions):
    rows, objects, differences = [], [], []
    for threads, directory in ((1, "one"), (4, "four")):
        for kernel in ("original", "retained"):
            paths = [root / directory / f"{version}-{kernel}" for version in versions]
            pair = [read_profile(path, version, kernel, threads) for path, version in zip(paths, versions)]
            (before, first_report, first_text), (after, last_report, last_text) = pair
            for field in ("binary_sha256", "architecture", "event"):
                if before[field] != after[field]:
                    raise ValueError(f"Different {field} for {kernel}, {threads} threads")
            rows.extend([before, after])
            diff = list(unified_diff(first_text.splitlines(), last_text.splitlines(),
                                     fromfile=versions[0], tofile=versions[1], lineterm=""))
            differences.append({"kernel": kernel, "threads": threads, "compared_versions": versions,
                                "same_hot_loop_id": before["loop"] == after["loop"],
                                "same_full_cqa_report": first_report == last_report,
                                "initial_advice_titles_by_path": titles(first_report), "full_report_diff": diff})
            (old, old_samples), (new, new_samples) = map(object_inventory, paths)
            for key in sorted(old.keys() | new.keys()):
                old_asm, old_text = old.get(key, ("", ""))
                new_asm, new_text = new.get(key, ("", ""))
                if not old_text or not new_text:
                    status = "previous_only" if old_text else "current_only"
                    report = json.loads(old_text or new_text)
                    differences[-1].setdefault("unpaired_reports", {})[key] = {
                        "version": versions[0] if old_text else versions[1],
                        "advice_titles_by_path": titles(report),
                    }
                elif old_asm != new_asm:
                    status = "different_instructions"
                else:
                    status = "identical" if old_text == new_text else "changed_report"
                objects.append({"kernel": kernel, "threads": threads, "object": key,
                                "type": "function" if key.startswith("fct_") else "loop", "status": status,
                                "previous_version": versions[0], "current_version": versions[1],
                                "previous_inclusive_sample_pct": old_samples.get(key, ""),
                                "current_inclusive_sample_pct": new_samples.get(key, ""),
                                "previous_assembly_sha256": old_asm, "current_assembly_sha256": new_asm,
                                "previous_report_sha256": sha256(old_text.encode()).hexdigest() if old_text else "",
                                "current_report_sha256": sha256(new_text.encode()).hexdigest() if new_text else ""})
                if status == "changed_report" and key != before["loop"]:
                    diff = unified_diff(old_text.splitlines(), new_text.splitlines(),
                                        fromfile=versions[0], tofile=versions[1], lineterm="")
                    differences[-1].setdefault("other_changed_reports", {})[key] = list(diff)
    return rows, objects, differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="directory containing one/ and four/ extracted experiments")
    parser.add_argument("output", type=Path, help="destination directory for diagnostics CSV files and JSON diffs")
    parser.add_argument("--versions", nargs=2, default=["2025.1.4", "2026.1.1"])
    args = parser.parse_args()
    rows, objects, differences = compare(args.root, args.versions)
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, values in (("diagnostics.csv", rows), ("diagnostics-objects.csv", objects)):
        with (args.output / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(values)
    (args.output / "diagnostics.json").write_text(json.dumps(differences, indent=2, ensure_ascii=False) + "\n")
    for case in differences:
        if not case["same_full_cqa_report"] or not case["same_hot_loop_id"]:
            print(f"Changed CQA report or hot loop: {case['kernel']} at {case['threads']} threads (diagnostics.json)")
    for obj in objects:
        if obj["status"] != "identical":
            print(f"{obj['object']} at {obj['threads']} threads: {obj['status']} (diagnostics-objects.csv)")


if __name__ == "__main__":
    main()
