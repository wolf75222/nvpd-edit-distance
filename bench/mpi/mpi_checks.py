#!/usr/bin/env python3
"""Compare MPI command output with compiled scalar and edlib oracles."""
import argparse
import json
import os
import random
import shlex
import subprocess
import tempfile
from pathlib import Path


def invoke(cmd):
    run = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
    if run.returncode:
        raise RuntimeError(f"{cmd} exited {run.returncode}\n{run.stderr[-1500:]}")
    return run.stdout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", required=True)
    p.add_argument("--oracle", required=True)
    p.add_argument("--launcher", required=True, help="e.g. 'mpiexec -n 4'")
    a = p.parse_args()
    launch = shlex.split(a.launcher) + [a.binary]
    rng = random.Random(87)
    pairs = [("", ""), ("", "AC"), ("A", "T"), ("TT", "R"), ("NNNNN", "G"),
             ("CCCCAAA", "AGAAAAACCAAAAAA"), ("NRGACNT", "RAAAGACNNNAACTTAA"),
             ("R" + "A"*11 + "C", "C" + "A"*33 + "G"*51 + "R"),
             ("RYNNR", "R" + "A"*511 + "Y"),
             ("CGT" + "A"*24, "A"*110 + "CGT" + "A"*27),
             ("RYSWKMBDHVN"*3, "NNN" + "RYSWKMBDHVN"*5 + "W")]
    pairs.extend(("".join(rng.choices("ACGTRYN", k=m)), "".join(rng.choices("ACGTRY", k=n)))
                 for m, n in [(7,13),(19,59),(33,17),(68,123)])
    with tempfile.TemporaryDirectory(dir=os.environ.get("CHPS_MPI_DATA_TMP")) as tmp:
        fasta = Path(tmp) / "input.fa"
        for i, (x,y) in enumerate(pairs):
            fasta.write_text(f">a\n{x}\n>b\n{y}\n")
            args = ["pair","--input",str(fasta)]
            native = [json.loads(invoke([a.oracle,*args,"--algorithm",algo]))["distance"]
                      for algo in ("seq","edlib")]
            height = (1, 3, 8, 32, 2**31-1)[i % 5]
            for extra in ([], ["--pipeline-rows", str(height)]):
                actual = json.loads(invoke([*launch, *args, "--threads", str(1+i%2), *extra]))
                if native != [actual["distance"]]*2:
                    raise ValueError(f"pair {i}, {extra}: scalar/edlib={native} MPI={actual}")
                method = "nvpd-mpi-pipeline" if extra else "nvpd-mpi"
                if actual["algorithm"] != method:
                    raise ValueError(f"wrong algorithm label: {actual}")
        for seqs in [["A"], ["", "ANC"], ["GAT", "TT", "", "ATC", "RNNN", "GAT"]]:
            fasta.write_text("".join(f">r{i}\n{s}\n" for i,s in enumerate(seqs)))
            args = ["matrix","--input",str(fasta)]
            expected = invoke([a.oracle,*args,"--algorithm","edlib"])
            omp_matrix = invoke([a.oracle, *args, "--algorithm", "nvpd-omp", "--threads", "4"])
            if expected != omp_matrix:
                raise ValueError("OpenMP matrix disagrees with edlib")
            actual = invoke([*launch,*args,"--threads","2"])
            if expected != actual:
                raise ValueError("MPI matrix disagrees with edlib")
    print(f"{len(pairs)} pairs in both MPI modes matched scalar and edlib, 3 MPI and 3 OpenMP matrices against edlib")


if __name__ == "__main__":
    main()
