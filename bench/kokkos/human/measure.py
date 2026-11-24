"""Check the true biological scores, then time independent warmed processes in matched blocks."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def references(binary, inputs, out):
    data = list(csv.DictReader((inputs / "inputs.csv").open()))
    for r in data:
        path = inputs / (r["case"] + ".fa")
        if hashlib.sha256(path.read_bytes()).hexdigest() != r["input_sha256"]:
            raise ValueError("FASTA digest mismatch")
        scores = []
        for alg in ("edlib", "wfa2"):
            obs = json.loads(
                subprocess.check_output(
                    [binary, "pair", "--algorithm", alg, "--input", path, "--threads", "1"]
                )
            )
            if int(r["m"]) != obs["m"] or int(r["n"]) != obs["n"]:
                raise ValueError("sequence length mismatch")
            scores.append(obs["distance"])
        if scores[0] != scores[1]:
            raise ValueError("independent library disagreement")
        r.update(distance=scores[0], normalized_distance=scores[0] / max(int(r["m"]), int(r["n"])))
    write(out / "references.csv", data)


def measure(target, binary, inputs, out, lane):
    workers = (1, 4, 8)[lane] if target == "epyc" else 1
    mask = sorted(os.sched_getaffinity(0))[:workers]
    if len(mask) != workers:
        raise ValueError("insufficient CPU affinity")
    ranks = ";".join(map(str, mask))
    host = socket.gethostname().split(".")[0]
    binary_sha = hashlib.sha256(binary.read_bytes()).hexdigest()
    fields = (
        "target job host case region algorithm execution_space workers cpu_ids block order m n distance "
        "seconds cells gcups peak_rss_kib host_workspace_array_bytes device_bytes input_sha256 binary_sha256"
    ).split()
    with (out / f"lane-{lane}.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for r in csv.DictReader((out / "references.csv").open()):
            m, n = int(r["m"]), int(r["n"])
            if target == "epyc" and m > 65536:
                continue
            path = inputs / (r["case"] + ".fa")
            if hashlib.sha256(path.read_bytes()).hexdigest() != r["input_sha256"]:
                raise ValueError("input changed after oracle")
            algs = ["native", "kokkos-rows", "kokkos-opt"]
            print(target, workers, r["case"], flush=True)
            for block in range(1, 6):
                k = (block - 1) % 3
                order = algs if block <= 3 else algs[::-1]
                for pos, alg in enumerate(order[k:] + order[:k], 1):
                    argv = [
                        "taskset",
                        "-c",
                        ",".join(map(str, mask)),
                        binary,
                        alg,
                        path,
                        str(workers),
                        r["distance"],
                    ]
                    line = subprocess.check_output(
                        argv, text=True, env=dict(os.environ, OMP_NUM_THREADS=str(workers))
                    ).strip()
                    space, score, sec, rss, dev = line.split(",")
                    if score != r["distance"] or not float(sec) > 0 or int(rss) <= 0:
                        raise ValueError(f"invalid score/time/RSS {line}")
                    # All manifest inputs contain the four canonical bases. The rolling
                    # table has two rows. Kokkos also makes a separate index view.
                    width = 4 * (max(m, n) + 1)
                    host_arrays = width * (4 if target == "gh200" else (6 if alg == "native" else 10))
                    writer.writerow(
                        dict(
                            target=target,
                            job=os.environ["SLURM_JOB_ID"],
                            host=host,
                            case=r["case"],
                            region=r["region"],
                            algorithm=alg,
                            execution_space=space,
                            workers=workers,
                            cpu_ids=ranks,
                            block=block,
                            order=pos,
                            m=m,
                            n=n,
                            distance=score,
                            seconds=sec,
                            cells=m * n,
                            gcups=m * n / float(sec) / 1e9,
                            peak_rss_kib=rss,
                            host_workspace_array_bytes=host_arrays,
                            device_bytes=dev,
                            input_sha256=r["input_sha256"],
                            binary_sha256=binary_sha,
                        )
                    )
                    stream.flush()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["reference", "epyc", "gh200"])
    p.add_argument("binary", type=Path)
    p.add_argument("inputs", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    if args.mode == "reference":
        references(args.binary, args.inputs, args.out)
    else:
        measure(args.mode, args.binary, args.inputs, args.out, int(os.environ.get("SLURM_LOCALID", 0)))


if __name__ == "__main__":
    main()
