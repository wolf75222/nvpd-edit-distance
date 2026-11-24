"""Measure time inside Allgather on the rank that determines the pair runtime."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import random
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mpi_benchmark import generated, launch, run, sha, sizes


def measure(binary, oracle, library, directory, output):
    pairs = {}
    for name, mode, path in generated(directory):
        if mode == "pair":
            score = json.loads(run([oracle, "pair", "--input", path, "--algorithm", "edlib"]).stdout)["distance"]
            pairs[name] = (path, score, sha(path.read_bytes()), min(sizes(path)))
    rng = random.Random(5634)
    fields = "job case nodes instrumented block distance input_sha256 physical_cores seconds calls critical_rank collective_seconds".split()
    with output.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for block in range(1, 8):
            configs = [(name, nodes, active) for name in pairs for nodes in (1, 2) for active in (False, True)]
            rng.shuffle(configs)
            for name, nodes, active in configs:
                path, score, digest, rows = pairs[name]
                cmd = launch(nodes, 24, 8)
                if active:
                    cmd += ["env", f"LD_PRELOAD={library}"]
                proc = run([*cmd, binary, "pair", "--input", path, "--threads", "8", "--warmups", "2", "--trials", "1"])
                logs = [json.loads(line) for line in proc.stderr.splitlines() if line.startswith("{")]
                [timing] = [row for row in logs if "trial" in row]
                measured = [row for row in logs if "collective" in row]
                if (timing["sampled_physical_cores"] != nodes*192 or json.loads(proc.stdout)["distance"] != score
                        or len(measured) != (3 if active else 0)):
                    raise ValueError(f"unexpected placement, distance or capture: {name} {nodes}")
                row = dict(job=os.environ["SLURM_JOB_ID"], case=name, nodes=nodes, instrumented=int(active),
                           block=block, distance=score, input_sha256=digest, physical_cores=nodes*192, seconds=timing["seconds"])
                if active:
                    timed = measured[-1]
                    if timed["calls"] != rows or not math.isclose(timed["total"], timing["seconds"], rel_tol=1e-8):
                        raise ValueError("unexpected rank or collective count")
                    row.update(calls=rows, critical_rank=timed["critical_rank"], collective_seconds=timed["collective"])
                writer.writerow(row)
                stream.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("binary", "oracle", "library", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(dir=os.environ.get("SLURM_SUBMIT_DIR", ".")) as work:
        measure(args.binary.resolve(), args.oracle.resolve(), args.library.resolve(), Path(work), args.output)
