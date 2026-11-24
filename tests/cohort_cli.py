"""Check literal-symbol, duplicate, and word-boundary scores at public CLI entry points."""
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def run(binary, *args):
    return subprocess.check_output([binary, *map(str, args)], text=True)


def main(binary, mpi=None):
    rng = random.Random(257)
    dna = "".join(rng.choices("ACGT", k=257))
    sequences = [dna[:n] for n in (1,7,8,31,32,33,63,64,65,128,257)]
    sequences += [dna, dna[1:], dna[:32]+"N"+dna[33:], "RYN", "N", "A"*192, "C"*192]
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp)/"input.fa"
        for panel in (sequences, [dna,dna,dna], ["N","A","N","R","Y", "NNNNNNR", "ANGTANGTANGT", "ARY", "RRG"]):
            path.write_text("".join(f">id{i}\n{s}\n" for i,s in enumerate(panel)))
            args = ["matrix", "--input", path, "--threads","2"]
            expected = run(binary,*args,"--algorithm","edlib")
            for kernel in ("packed-wfa","nvpd-sparse","nvpd-split","wfa2"):
                for options in ([],["--unique"]):
                    if run(binary,*args,"--algorithm",kernel,*options) != expected:
                        raise ValueError(f"{kernel} failed to restore matrix labels or literal scores")
            if mpi:
                rows = [line.split("\t")[1:] for line in expected.splitlines()[1:]]
                certificates = sum(
                    max(len(a), len(b)) // 3 >= min(len(a), len(b))
                    and int(rows[i][j]) == abs(len(a) - len(b))
                    for i, a in enumerate(panel) for j, b in enumerate(panel[:i])
                )
                for kernel in ("nvpd-indexed","nvpd-omp","edlib","wfa2"):
                    observed = subprocess.run([mpi,*map(str,args),"--algorithm",kernel,"--length-bound","1",
                                               "--grain","1000"],
                                              capture_output=True,text=True,check=True)
                    if observed.stdout != expected: raise ValueError(f"MPI {kernel} disagrees with edlib")
                    diagnostics = [json.loads(line) for line in observed.stderr.splitlines() if line.startswith("{")]
                    ranks = [r for r in diagnostics if "rank" in r]
                    timed = [r for r in diagnostics if "mode" in r]
                    if len(ranks) != 1 or ranks[0]["rank"] != 0 or ranks[0]["certified_pairs"] != certificates:
                        raise ValueError(f"MPI {kernel} counted the wrong number of subsequences")
                    if len(timed) != 1 or timed[0]["assigned_workers"] != 1:
                        raise ValueError("one matrix chunk must have one active worker")
        path.write_text(f">a\n{dna}\n>b\nN{dna[3:]}T\n")
        for kernel in ("packed-wfa","nvpd-sparse"):
            args = ["pair","--input",path,"--algorithm",kernel]
            score = json.loads(run(binary,*args))["distance"]
            path_score = json.loads(run(binary,*args,"--trace"))
            if path_score["distance"] != score or not path_score["cigar"]:
                raise ValueError(f"{kernel} trace and score differ")


if __name__ == "__main__":
    main(*sys.argv[1:])
