"""Exercise matrix placement on small, empty and duplicate records without biological data."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

POLICIES = ("static","cyclic","hilbert","steal","hilbert-weighted")


def run(cmd):
    return subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=45)


def check(solver, oracle, workers="openmp"):
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/"input.fa"
        for seqs in ([""],["AR","","AR","NNRY","S"],["AA","A"],["","","",""]):
            path.write_text("".join(f">id{i}\n{x}\n" for i,x in enumerate(seqs)))
            expected=run([oracle,"matrix","--algorithm","edlib","--input",path]).stdout
            for schedule in POLICIES:
                args=["mpirun","--oversubscribe","-n","3",solver,"matrix","--input",path,
                      "--schedule",schedule,"--workers",workers,"--threads","2","--grain","4","--warmups","1"]
                result=run(args)
                assert result.stdout==expected
                ranks=[json.loads(x) for x in result.stderr.splitlines() if x.startswith('{"rank"')]
                assert len(ranks)==3
                assert sum(r["pairs"] for r in ranks)==len(seqs)*(len(seqs)-1)//2
                assert sum(r["cells"] for r in ranks)==sum(len(a)*len(b) for i,a in enumerate(seqs) for b in seqs[i+1:])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("solver",type=Path)
    p.add_argument("oracle",type=Path)
    p.add_argument("--workers",choices=("openmp","kokkos"),default="openmp")
    opts=p.parse_args()
    check(opts.solver.resolve(),opts.oracle.resolve(),opts.workers)
