"""Time matched chromosome windows and compare every score to global edlib."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mpi_benchmark import launch, run, run_rank_zero


def digest(data):
    return hashlib.sha256(data).hexdigest()


def matched_scores(path, oracle, temp):
    lines = path.read_bytes().splitlines()
    if len(lines) % 4:
        raise ValueError("window file must contain adjacent pairs on single sequence lines")
    expected, cells = [], 0
    for i in range(0, len(lines), 4):
        a, x, b, y = lines[i:i+4]
        temp.write_bytes(b"\n".join((a, x, b, y))+b"\n")
        score = json.loads(run([oracle,"pair","--input",temp,"--algorithm","edlib"]).stdout)["distance"]
        expected.append((a[1:].decode(),b[1:].decode(),score))
        cells += len(x)*len(y)
    return expected, cells


def describe(binary, oracle, inputs, out, work):
    cases = {}
    with (out/"scores.tsv").open("x",newline="") as stream:
        scores = csv.writer(stream,delimiter="\t",lineterminator="\n")
        scores.writerow(("id_a","id_b","edlib"))
        for key in ("matched-768","matched-384","pair-64k","pair-128k"):
            f = inputs/(key+".fa")
            if key=="matched-384":
                expected = cases["matched-768"][1][:384]
                big = (inputs/"matched-768.fa").read_bytes()
                if not big.startswith(f.read_bytes()) or len(f.read_bytes().splitlines())!=4*384:
                    raise ValueError("the 384-pair cohort must be a prefix of the 768-pair cohort")
                data=f.read_bytes().splitlines()
                cells=sum(len(data[i+1])*len(data[i+3]) for i in range(0,len(data),4))
            else:
                expected, cells = matched_scores(f,oracle,work/"single.fa")
                scores.writerows(expected)
            cases[key] = f, expected, cells, digest(f.read_bytes())

    configs=[("matched-384",1,0),("matched-768",1,0),("matched-768",2,0)]
    configs.extend((k,n,h) for k in ("pair-64k","pair-128k") for n in (1,2) for h in (0,32,128))
    rng = random.Random(18533)
    fields="job block order case nodes ranks threads physical_cores height pairs cells seconds gcups max_rank_rss_kib input_sha256 scores_sha256".split()
    with (out/"trials.csv").open("x",newline="") as stream:
        wr=csv.DictWriter(stream,fields)
        wr.writeheader()
        for block in range(1,6):
            rng.shuffle(configs)
            for order,(key,nodes,h) in enumerate(configs,1):
                f, expected, cells, fingerprint = cases[key]
                args=[binary,"pairs" if key.startswith("matched") else "pair","--input",f,"--threads","8",
                      "--warmups","1","--trials","1","--report-rss","1"]
                if h: args.extend(["--pipeline-rows",str(h)])
                result = run_rank_zero(launch(nodes,24,8),args,work)
                if len(expected)>1:
                    observed=[(r["id_a"],r["id_b"],int(r["distance"])) for r in csv.DictReader(io.StringIO(result.stdout),delimiter="\t")]
                else:
                    observed=[(*expected[0][:2],json.loads(result.stdout)["distance"])]
                if observed!=expected: raise ValueError(f"edlib discrepancy in {key}")
                [t]=[json.loads(x) for x in result.stderr.splitlines() if x.startswith('{"mode"')]
                if (t["assigned_workers"],t["sampled_physical_cores"],t["ranks"],t["pipeline_rows"],t["trial"])!=(192*nodes,192*nodes,24*nodes,h,1):
                    raise ValueError(f"invalid worker count or method for {key}: {t}")
                if t["seconds"]<=0 or t["max_rank_rss_kib"]<=0: raise ValueError("missing time or RSS")
                row=dict(job=os.environ["SLURM_JOB_ID"],block=block,order=order,case=key,nodes=nodes,ranks=24*nodes,
                         threads=8,physical_cores=192*nodes,height=h,pairs=len(expected),cells=cells,seconds=t["seconds"],
                         gcups=cells/t["seconds"]/1e9,max_rank_rss_kib=t["max_rank_rss_kib"],input_sha256=fingerprint,
                         scores_sha256=digest("\n".join(str(s[2]) for s in expected).encode()))
                wr.writerow(row)
                stream.flush()
            print(f"completed randomized block {block}",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("binary","oracle","inputs","out","work"):p.add_argument(name,type=Path)
    a=p.parse_args()
    describe(*(getattr(a,k).resolve() for k in ("binary","oracle","inputs","out","work")))
