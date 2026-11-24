#!/usr/bin/env python3
"""Weak scaling of exact, distributed edit distance, with a fixed DP workload per core."""
import argparse
import csv
import json
import os
import random
import tempfile
from pathlib import Path
from statistics import median

from mpi_benchmark import launch, run, sha

CORES = 192
TRIALS = 5


def inputs():
    rng = random.Random(192384)

    def dna(n):
        return "".join(rng.choices("ACGT", k=n))

    def mutate(sequence):
        return "".join(rng.choice("ACGT".replace(base, "")) if rng.random() < 0.02 else base
                       for base in sequence)

    rows = dna(1024)
    columns = mutate(rows) + dna(2 * CORES * 4096 - len(rows))
    for width in (256, 4096):
        yield f"pair_{width}", "pair", [[rows, columns[:CORES * nodes * width]] for nodes in (1, 2)]
    ancestor = dna(8192)
    cohort = [mutate(ancestor) for _ in range(120)]
    yield "matrix", "matrix", [cohort[:85], cohort]


def write_table(path, rows):
    with path.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def benchmark(args, work):
    trials, summaries = [], []
    for case_num, (case, mode, layouts) in enumerate(inputs()):
        work_per_core, times = {}, {}
        for nodes in ((1, 2) if case_num % 2 == 0 else (2, 1)):
            seqs = layouts[nodes - 1]
            path = work / "input.fa"
            path.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(seqs)))
            fasta_sha = sha(path.read_bytes())
            cells = sum(len(x) * len(y) for i, x in enumerate(seqs) for y in seqs[:i])
            cores = nodes * CORES
            work_per_core[nodes] = cells / cores
            options = [mode, "--input", str(path), "--threads", "8"]
            oracle = run([str(args.oracle), *options, "--algorithm", "edlib"]).stdout
            proc = run([*launch(nodes, 24, 8), str(args.binary), *options, "--warmups", "1",
                        "--trials", str(TRIALS)])
            if mode == "pair":
                got = json.loads(proc.stdout)["distance"]
                equal = got == json.loads(oracle)["distance"]
            else:
                got = sha(proc.stdout.encode())
                equal = proc.stdout == oracle
            stats = [json.loads(x) for x in proc.stderr.splitlines() if x.startswith('{"mode":')]
            if not equal or [r["trial"] for r in stats] != list(range(1, TRIALS + 1)):
                raise RuntimeError(f"oracle or trial mismatch for {case} on {nodes} nodes")
            for r in stats:
                if (r["ranks"], r["threads_per_rank"], r["assigned_workers"],
                    r["sampled_physical_cores"]) != (24 * nodes, 8, cores, cores):
                    raise RuntimeError(f"unexpected physical placement for {case}: {r}")
                trials.append(dict(job=os.environ["SLURM_JOB_ID"], case=case, nodes=nodes, **r,
                                   reservation_cores=384, records=len(seqs), first_length=len(seqs[0]),
                                   last_length=len(seqs[-1]), pairs=len(seqs)*(len(seqs)-1)//2,
                                   cells=cells, cells_per_core=work_per_core[nodes],
                                   gcups=cells/r["seconds"]/1e9, input_sha256=fasta_sha,
                                   oracle_algorithm="edlib", distance_or_matrix_sha256=got))
            times[nodes] = median(r["seconds"] for r in stats)
            print(f"{case} {nodes} node(s), {times[nodes]:.6f} s", flush=True)
        if work_per_core[1] != work_per_core[2]:
            raise RuntimeError("cells/core differ across node counts")
        efficiency = times[1] / times[2]
        summaries.append(dict(case=case, median_seconds_1_node=times[1], median_seconds_2_nodes=times[2],
                              weak_efficiency=efficiency, throughput_ratio=2*efficiency,
                              cells_per_core=work_per_core[1]))
    write_table(args.out / "trials.csv", trials)
    write_table(args.out / "summary.csv", summaries)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--oracle", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=os.environ.get("CHPS_MPI_DATA_TMP")) as tmp:
        benchmark(args, Path(tmp))
