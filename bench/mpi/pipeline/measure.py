"""Paired row-wise vs nearest-neighbor NvPD, with seven randomized process blocks."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import random
import statistics
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mpi_benchmark import generated, launch, run, run_rank_zero, sha, sizes

def compare(binary, oracle, probe, directory, output):
    cases = {}
    for name, mode, path in generated(directory):
        if mode == "pair":
            result = run([oracle, "pair", "--input", path, "--algorithm", "edlib"])
            m, n = sorted(sizes(path))
            cases[name] = (path, json.loads(result.stdout)["distance"], sha(path.read_bytes()), m, n)
    previous_sha = {"square": "03ea4833159dd1374c752496592e5a3008658f7a42f2d09f08af97ffa80b10d2",
                    "wide": "8790d9c24414cfc977cd1107c80ccd647ae25e02aad6e3fb90fd340f9f47b0be"}
    if {name: data[2] for name, data in cases.items()} != previous_sha:
        raise ValueError("inputs differ from the historical MPI campaign")
    rng = random.Random(41721)
    configs = [(c, nodes, h) for c in cases for nodes in (1, 2) for h in (0, 8, 32)]
    records = []
    with (output / "trials.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, "case nodes height block order ranks threads physical_cores m n distance sha256 seconds".split())
        writer.writeheader()
        with (output / "allgather.csv").open("x", newline="") as events:
            counts = csv.writer(events)
            counts.writerow(("case", "nodes", "height", "allgathers_per_rank"))
            for block in range(8):
                rng.shuffle(configs)
                for order, (name, nodes, height) in enumerate(configs, 1):
                    path, exact, digest, m, n = cases[name]
                    cmd = launch(nodes, 24, 8)
                    if block == 0:
                        cmd += ["env", f"LD_PRELOAD={probe}"]
                    arguments = [binary, "pair", "--input", path, "--threads", "8",
                            "--warmups", "2" if block else "0", "--trials", "1"]
                    if height:
                        arguments += ["--pipeline-rows", str(height)]
                    proc = run_rank_zero(cmd, arguments, directory)
                    result = json.loads(proc.stdout)
                    logs = [json.loads(line) for line in proc.stderr.splitlines() if line.startswith("{")]
                    [t] = [entry for entry in logs if "trial" in entry]
                    method = "nvpd-mpi-pipeline" if height else "nvpd-mpi"
                    if result != {"algorithm": method, "distance": exact}:
                        raise ValueError((name, height, result))
                    if (t["pipeline_rows"], t["ranks"], t["threads_per_rank"], t["sampled_physical_cores"],
                            t["assigned_workers"]) != (height, 24*nodes, 8, 192*nodes, 192*nodes):
                        raise ValueError(f"incorrect placement or mode {t}")
                    if block == 0:
                        [p] = [entry for entry in logs if "collective" in entry]
                        if p["calls"] != (0 if height else m):
                            raise ValueError(f"unexpected Allgather count: {p}")
                        counts.writerow([name, nodes, height, p["calls"]])
                        events.flush()
                        continue
                    r = dict(case=name, nodes=nodes, height=height, block=block, order=order,
                             ranks=24*nodes, threads=8, physical_cores=192*nodes, m=m, n=n,
                             distance=exact, sha256=digest, seconds=t["seconds"])
                    records.append(r)
                    writer.writerow(r)
                    stream.flush()
    with (output / "summary.csv").open("x", newline="") as f:
        w = csv.writer(f)
        w.writerow("case nodes height units metric mean sd low95 high95".split())
        times = {(r["case"], r["nodes"], r["height"], r["block"]): r["seconds"] for r in records}
        for name, nodes, h in sorted(configs):
            keys = range(1, 8)
            wall = [times[name, nodes, h, b] for b in keys]
            gain = [times[name, nodes, 0, b] / times[name, nodes, h, b] for b in keys]
            scale = [times[name, 1, h, b] / times[name, nodes, h, b] for b in keys]
            for metric, group in (("seconds", wall), ("speedup_vs_allgather", gain), ("speedup_vs_1node", scale)):
                avg, sd = statistics.mean(group), statistics.stdev(group)
                half = 2.446911851 * sd / math.sqrt(len(group))
                w.writerow((name, nodes, h, len(group), metric, avg, sd, avg-half, avg+half))

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("binary", "oracle", "probe", "output"):
        p.add_argument(key, type=Path)
    args = p.parse_args()
    with tempfile.TemporaryDirectory(dir=os.environ.get("SLURM_SUBMIT_DIR", ".")) as temp:
        compare(args.binary.resolve(), args.oracle.resolve(), args.probe.resolve(), Path(temp), args.output)
