#!/usr/bin/env python3
"""Exercise the public Kokkos CLI and record same-binary, warmed A/B trials."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

def invoke(cmd, timeout=120):
    return subprocess.check_output([str(x) for x in cmd], text=True, timeout=timeout)


def check_cli(binary, work, threads):
    pair = work / "edge.fa"
    edges = [("", ""), ("", "N"), ("AA", ""), ("ACGTRY", "NGY"),
             ("T" * 7 + "R", "A" * 260 + "R"), ("TGAC", "TGAC"), ("GA", "A" * 259)]
    for a, b in edges:
        pair.write_text(f">a\n{a}\n>b\n{b}\n")
        args = [binary, "pair", "--input", pair, "--threads", threads, "--algorithm"]
        expected = json.loads(invoke([*args, "seq", "--trace"]))
        cost = json.loads(invoke([*args, "nvpd-kokkos"]))
        trace = json.loads(invoke([*args, "nvpd-kokkos", "--trace"]))
        if expected["distance"] != cost["distance"] or expected["distance"] != trace["distance"]:
            raise ValueError("CLI cost/trace disagrees with sequential distance")
        if expected["cigar"] != trace["cigar"]:
            raise ValueError("traceback differs from deterministic serial traceback")
    seqs = ["", "AN", "T", "AAAA", "CTRG"]
    pair.write_text("".join(f">tip{i}\n{a}\n" for i, a in enumerate(seqs)))
    cmd = [binary, "matrix", "--input", pair, "--threads", threads, "--algorithm"]
    if invoke([*cmd, "nvpd-kokkos"]) != invoke([*cmd, "edlib"]):
        raise ValueError("Kokkos pairwise matrix differs from edlib")
    return {"cost_and_trace_pairs": len(edges), "matrix_tips": len(seqs),
            "unique_matrix_pairs": len(seqs) * (len(seqs) - 1) // 2}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("target", choices=("cpu", "gpu"))
    p.add_argument("build", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("threads", type=int)
    a = p.parse_args()
    binary = a.build / "editdist"
    runner = a.build / "measure_kokkos"
    checks = check_cli(binary, a.work, a.threads)
    (a.out / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    data = a.work / "input"
    generator = Path(__file__).resolve().parents[1] / "scripts/synthetic.py"
    invoke([sys.executable, generator, "--out", data, "--lengths", 5000, 20000, "--rates", 0.05, "--pairs", 1])
    with (data / "manifest.tsv").open() as f:
        pairs = list(csv.DictReader(f, delimiter="\t"))
    native = "nvpd-omp" if a.target == "cpu" else "nvpd-cuda"
    fields = ["target", "execution_space", "algorithm", "threads", "m", "n", "input_sha256", "trial",
              "seconds", "peak_rss_kib", "device_bytes", "distance", "reference_distance", "job"]
    workers = sorted({1, a.threads}) if a.target == "cpu" else [a.threads]
    job = os.environ.get("SLURM_JOB_ID", "local")
    with (a.out / "timings.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for pair in pairs:
            fasta = data / pair["path"]
            digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
            edlib = json.loads(invoke([binary, "pair", "--algorithm", "edlib", "--input", fasta]))["distance"]
            for threads in workers:
                for trial in range(1, 4):
                    algs = [native, "nvpd-kokkos"] if trial % 2 else ["nvpd-kokkos", native]
                    for alg in algs:
                        raw = json.loads(invoke([runner, alg, fasta, threads]))
                        if raw["distance"] != edlib:
                            raise ValueError(f"{alg} disagrees with edlib at {pair['length']} bases")
                        if alg == "nvpd-kokkos" and raw["execution_space"].lower() != (
                                "openmp" if a.target == "cpu" else "cuda"):
                            raise ValueError("Kokkos executed in an unexpected space")
                        writer.writerow(dict(raw, target=a.target, algorithm=alg, threads=threads,
                                             input_sha256=digest, trial=trial, reference_distance=edlib, job=job))
                        stream.flush()
    human = os.environ.get("CHPS_HUMAN_PAIR")
    if a.target == "gpu" and human:
        expected = int(os.environ["CHPS_HUMAN_EXPECTED"])
        try:
            res = json.loads(invoke([binary, "pair", "--input", human, "--algorithm", "nvpd-kokkos",
                                     "--threads", a.threads], int(os.environ.get("CHPS_KOKKOS_TIMEOUT", 180))))
            res.update(reference_distance=expected, status="ok" if res["distance"] == expected else "mismatch",
                       input_sha256=hashlib.sha256(Path(human).read_bytes()).hexdigest(), job=job)
        except subprocess.TimeoutExpired:
            res = {"status": "timeout", "reference_distance": expected, "job": job}
        (a.out / "human.json").write_text(json.dumps(res, indent=2) + "\n")
        if res["status"] == "mismatch": raise ValueError("human distance mismatch")


if __name__ == "__main__":
    main()
