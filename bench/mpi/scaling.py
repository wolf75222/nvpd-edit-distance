"""Repeat strong and weak scaling in separately warmed MPI launches."""
import argparse
import csv
import json
import os
from pathlib import Path
import random
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mpi_benchmark import generated, launch, run, sha, sizes
from mpi_weak import inputs


def scaling(binary, oracle, cohort, directory, output):
    cases = {}
    for name, mode, path in generated(directory) + [("phix64", "matrix", cohort)]:
        for nodes in (1, 2):
            cases["strong", name, nodes] = (mode, path)
    for name, mode, layouts in inputs():
        for nodes, seqs in enumerate(layouts, 1):
            path = directory / f"{name}-{nodes}.fa"
            path.write_text("".join(f">r{i}\n{s}\n" for i,s in enumerate(seqs)))
            cases["weak", name, nodes] = (mode, path)
    metadata = {}
    for mode, path in set(cases.values()):
        lengths = sizes(path)
        pairs = len(lengths)*(len(lengths)-1)//2
        cells = (sum(lengths)**2-sum(x*x for x in lengths))//2
        ref = run([oracle,mode,"--input",path,"--algorithm","edlib","--threads","8"]).stdout
        score = str(json.loads(ref)["distance"]) if mode == "pair" else sha(ref.encode())
        metadata[path] = (score, pairs, cells, sha(path.read_bytes()))
    rng = random.Random(171)
    fields = "job study case nodes block ranks threads physical_cores pairs cells input_sha256 oracle seconds".split()
    with output.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for block in range(1, 8):
            keys = list(cases)
            rng.shuffle(keys)
            for study, name, nodes in keys:
                mode, path = cases[study,name,nodes]
                expected, pairs, cells, digest = metadata[path]
                # Capture the root streams at their source. The MPI launcher forwards
                # rank output asynchronously and may merge stdout and stderr.
                stdout, stderr = directory / "rank0.out", directory / "rank0.err"
                redirect = ('if [ "${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-}}" = 0 ]; then '
                            'exec >"$1" 2>"$2"; fi; shift 2; exec "$@"')
                run([*launch(nodes,24,8), "sh", "-c", redirect, "rank", stdout, stderr, binary,mode,"--input",path,"--threads","8",
                            "--warmups","2","--trials","1"])
                score = str(json.loads(stdout.read_text())["distance"]) if mode == "pair" else sha(stdout.read_bytes())
                [timing] = [json.loads(x) for x in stderr.read_text().splitlines() if x.startswith('{"mode":')]
                if (timing["sampled_physical_cores"], timing["assigned_workers"], score) != (nodes*192, nodes*192, expected):
                    raise ValueError(f"invalid score or worker placement on {name}/{nodes}")
                writer.writerow(dict(job=os.environ["SLURM_JOB_ID"], study=study, case=name, nodes=nodes, block=block,
                                     ranks=24*nodes, threads=8, physical_cores=192*nodes, pairs=pairs, cells=cells,
                                     input_sha256=digest, oracle=score, seconds=timing["seconds"]))
                stream.flush()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for field in ("binary", "oracle", "cohort", "output"):
        p.add_argument(field, type=Path)
    opts = p.parse_args()
    with tempfile.TemporaryDirectory(dir=os.environ.get("SLURM_SUBMIT_DIR", ".")) as work:
        scaling(opts.binary.resolve(), opts.oracle.resolve(), opts.cohort.resolve(), Path(work), opts.output)
