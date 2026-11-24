#!/usr/bin/env python3
"""Matched one- and two-node strong scaling for distributed NvPD."""
import argparse
import csv
import hashlib
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path
from statistics import median


def sha(data):
    return hashlib.sha256(data).hexdigest()


def generated(work):
    rng = random.Random(20261124)
    def dna(n):
        return "".join(rng.choices("ACGT", k=n))
    def mutate(s):
        return "".join(rng.choice("ACGT".replace(c,"")) if rng.random() < .02 else c for c in s)
    def fasta(name, records):
        path = work / (name+".fa")
        path.write_text("".join(f">r{i}\n{s}\n" for i,s in enumerate(records)))
        return path
    square = dna(20000)
    wide = dna(1024)
    anchor = dna(12000)
    return [("square", "pair", fasta("square",[square, mutate(square)+dna(37)])),
            ("wide", "pair", fasta("wide",[wide, mutate(wide)+dna(1000003-1024)])),
            ("synthetic64", "matrix", fasta("synthetic64",[mutate(anchor) for _ in range(64)]))]


def sizes(path):
    lengths = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            lengths.append(0)
        else:
            lengths[-1] += len(line.strip())
    return lengths


def launch(nodes, ranks_per_node, threads):
    if os.getenv("CHPS_MPI_LAUNCHER") == "mpirun":
        hosts = os.environ["CHPS_MPI_HOSTS"].split(",")[:nodes]
        return ["mpirun","--host",",".join(h+":192" for h in hosts),"-n",str(nodes*ranks_per_node),
                "--map-by",f"ppr:{ranks_per_node}:node:PE={threads}","--bind-to","core","--oversubscribe"]
    return ["srun","--exclusive","--exact","--mpi="+os.getenv("CHPS_SLURM_MPI","pmix"),
            "--hint=nomultithread","--cpu-bind=cores","--distribution=block:block",
            "-N",str(nodes),"-n",str(nodes*ranks_per_node),
            "--ntasks-per-node="+str(ranks_per_node),"-c",str(threads)]


def run(cmd):
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=420)
    if proc.returncode:
        raise RuntimeError(f"{cmd} exited {proc.returncode}\n{proc.stderr[-3500:]}")
    return proc


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--oracle", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--cohort", type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fields = ["job","case","mode","nodes","ranks","threads_per_rank","assigned_workers","sampled_physical_cores",
              "step_cores","reservation_cores","records","pairs","cells","trial","seconds",
              "gcups","input_sha256","oracle_algorithm","distance_or_matrix_sha256"]
    grouped = {}
    with tempfile.TemporaryDirectory(dir=os.environ.get("CHPS_MPI_DATA_TMP")) as tmp, (a.out/"trials.csv").open("x",newline="") as file:
        out = csv.DictWriter(file, fields)
        out.writeheader()
        cases = generated(Path(tmp))
        if a.cohort: cases.append(("phix64","matrix",a.cohort))
        # Equal configurations in one exclusive allocation. Alternate node count
        # order by dataset to reduce order bias between the two strong-scaling runs.
        for case_no,(case,mode,path) in enumerate(cases):
            dims = sizes(path)
            cells = sum(dims[i]*dims[j] for i in range(len(dims)) for j in range(i))
            oracle = run([str(a.oracle),mode,"--input",str(path),"--algorithm","edlib","--threads","8"]).stdout
            expected = json.loads(oracle)["distance"] if mode == "pair" else sha(oracle.encode())
            for nodes in ((1,2) if case_no%2==0 else (2,1)):
                ranks_per_node, threads = 24,8
                step = launch(nodes,ranks_per_node,threads)
                cmd = [*step,str(a.binary),mode,"--input",str(path),"--threads",str(threads),
                       "--warmups","1","--trials","3"]
                proc = run(cmd)
                got = json.loads(proc.stdout)["distance"] if mode == "pair" else sha(proc.stdout.encode())
                if got != expected: raise ValueError(f"{case} disagrees with edlib")
                stats = [json.loads(x) for x in proc.stderr.splitlines() if x.startswith('{"mode":')]
                if [r["trial"] for r in stats] != [1,2,3]: raise ValueError(f"missing MPI trials for {case}")
                for r in stats:
                    if (r["ranks"] != nodes*ranks_per_node or r["assigned_workers"] != 192*nodes
                            or r["sampled_physical_cores"] != nodes*192):
                        raise ValueError(f"worker team mismatch for {case}: {r}")
                    out.writerow(dict(r,job=os.environ["SLURM_JOB_ID"],case=case,nodes=nodes,
                                      step_cores=192*nodes,reservation_cores=384,records=len(dims),
                                      pairs=len(dims)*(len(dims)-1)//2,cells=cells,gcups=cells/r["seconds"]/1e9,
                                      input_sha256=sha(path.read_bytes()),oracle_algorithm="edlib",
                                      distance_or_matrix_sha256=got))
                file.flush()
                grouped[case,nodes] = median(r["seconds"] for r in stats)
                print(f"{case} {nodes} node(s), {grouped[case,nodes]:.6f} s",flush=True)
    with (a.out/"summary.csv").open("x",newline="") as f:
        out = csv.writer(f)
        out.writerow(["case","median_seconds_1_node","median_seconds_2_nodes","speedup_2_vs_1","efficiency_2_vs_1"])
        for case,_,_ in cases:
            one,two = grouped[case,1],grouped[case,2]
            out.writerow([case,one,two,one/two,one/(2*two)])


if __name__ == "__main__":
    main()
