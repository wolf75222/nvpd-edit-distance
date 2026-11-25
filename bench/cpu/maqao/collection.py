"""Export ONEVIEW collection checks from an archive unpacked next to its log."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys


def observation(profile):
    run = profile / "shared/run_0"
    reports = {}
    for level in (0, 1):
        summary = json.loads((run / f"summary_{level}.json").read_text())
        reports.update((item["name"], item) for item in summary["report"])
    log = Path(str(profile) + ".log").read_text()
    application_log = (profile / "logs/run_0/lprof.log").read_text()
    scores = re.findall(r"^(\d+),[\d.eE+-]+$", application_log, flags=re.M)
    if set(scores) != {"1479"}:
        raise ValueError(f"Unexpected scores in {profile}")
    warnings = [line.split(" Warning: ", 1)[1] for line in log.splitlines() if " Warning: " in line]
    return {
        "protocol": profile.parent.name, "profile": profile.name, "calls": len(scores),
        "score": scores[0], "duration_check": reports["lvl0_profiled_time"]["status"],
        "duration_message": reports["lvl0_profiled_time"]["text"],
        "topology_check": reports["lvl0_is_lstopo"]["status"],
        "debug_info_check": reports.get("lvl0_source_debug", {}).get("status", ""),
        "compiler_options": reports["lvl0_utilitary_options"]["text"],
        "mlock_buffer_warning": int(any("perf_event_mlock_kb" in w for w in warnings)),
        "kernel_sampling_warning": int(any("Kernel-space measurements not allowed" in w for w in warnings)),
        "level1_alerts": " | ".join(key for key, item in reports.items()
                                    if key.startswith("lvl1_") and item["status"] != "good"),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", nargs="+", type=Path)
    args = parser.parse_args()
    rows = [observation(p) for p in args.profiles]
    writer = csv.DictWriter(sys.stdout, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
