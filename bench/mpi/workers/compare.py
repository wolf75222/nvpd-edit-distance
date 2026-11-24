"""Compare the two host worker pools with the same rank assignment and NvPD kernel."""
import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "schedules"))
from measure import POLICIES, check_output, launch, run_rank_zero


def run(binary, fasta, references, out, workdir):
    lines = fasta.read_text().splitlines()
    ids = [s[1:] for s in lines if s.startswith(">")]
    with references.open() as stream:
        expected = {(r["id_a"],r["id_b"]):int(r["edlib"]) for r in csv.DictReader(stream,delimiter="\t")
                    if r["edlib"] == r["wfa2"]}
    if len(expected) != len(ids)*(len(ids)-1)//2:
        raise ValueError("incomplete reference matrix")
    configs = [(policy,nodes,backend) for policy in POLICIES for nodes in (1,2)
               for backend in ("kokkos","openmp")]
    with (out/"trials.csv").open("x") as stream:
        wr = csv.writer(stream)
        wr.writerow(("block","schedule","nodes","backend","physical_cores","seconds","peak_rank_rss_kib","pairs","cells"))
        for block in range(1,6):
            random.Random(19432 + block).shuffle(configs)
            for policy,nodes,backend in configs:
                command = [binary,"matrix","--input",fasta,"--schedule",policy,"--workers",backend,
                           "--grain","4","--threads","8","--warmups","1","--trials","1","--report-rss","1"]
                result = run_rank_zero(launch(nodes,24,8),command,workdir)
                check_output(result.stdout,ids,expected)
                diagnostics = [json.loads(line) for line in result.stderr.splitlines() if line.startswith("{")]
                timing, = [r for r in diagnostics if "mode" in r]
                ranks = [r for r in diagnostics if "rank" in r]
                pairs = sum(r["pairs"] for r in ranks)
                if (timing["worker_backend"] != backend or timing["sampled_physical_cores"] != 192*nodes
                        or pairs != len(expected) or len(ranks) != nodes*24):
                    raise ValueError(timing)
                wr.writerow((block,policy,nodes,backend,timing["sampled_physical_cores"],timing["seconds"],
                             timing["max_rank_rss_kib"],pairs,sum(r["cells"] for r in ranks)))
                stream.flush()


if __name__ == "__main__":
    run(*(Path(arg).resolve() for arg in sys.argv[1:]))
