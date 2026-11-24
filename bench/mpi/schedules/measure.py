"""Measure identical all-vs-all distances under five MPI scheduling policies."""
import argparse
import csv
import hashlib
import io
import itertools
import json
import os
from pathlib import Path
import random
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mpi_benchmark import launch, run_rank_zero
POLICIES = ("static", "cyclic", "hilbert", "steal", "hilbert-weighted")
TRIAL_FIELDS = ("job block order layout schedule nodes ranks threads cores grain pairs cells seconds gcups "
                "busy_fraction max_mean_rank_cells stolen_pairs cross_node_pairs max_rank_rss_kib "
                "input_sha256 canonical_scores_sha256").split()
WORK_FIELDS = "job block layout schedule nodes rank pairs cells stolen_pairs cross_node_pairs busy_seconds".split()


def manifest(path):
    with path.open() as f:
        records = list(csv.DictReader(f,delimiter="\t"))
    return {r["id"]:int(r["length"]) for r in records}


def check_output(text, ids, expected):
    table=csv.reader(io.StringIO(text),delimiter="\t")
    if next(table) != ["id",*ids]:
        raise ValueError("matrix header disagrees with input order")
    for label in ids:
        row=next(table)
        if row[0]!=label or len(row)!=len(ids)+1:
            raise ValueError("invalid matrix layout")
        for other,score in zip(ids,row[1:]):
            key=tuple(sorted((label,other)))
            if int(score)!=(0 if label==other else expected[key]):
                raise ValueError(f"incorrect matrix distance: {key}")
    if next(table,None) is not None:
        raise ValueError("unexpected matrix suffix")


def measure(binary,inputs,out,work,windows):
    lengths=manifest(windows)
    with (out/"scores.tsv").open() as f:
        refs=list(csv.DictReader(f,delimiter="\t"))
    expected={(r["id_a"],r["id_b"]):int(r["edlib"]) for r in refs if r["edlib"]==r["wfa2"]}
    if set(expected)!=set(itertools.combinations(sorted(lengths),2)):
        raise ValueError("incomplete independent score map")
    canonical=hashlib.sha256("".join(f"{a}\t{b}\t{d}\n" for (a,b),d in sorted(expected.items())).encode()).hexdigest()
    cells=sum(lengths[a]*lengths[b] for a,b in expected)
    layouts={}
    for name in ("locus","length"):
        blob=(inputs/(name+".fa")).read_bytes()
        ids=[x[1:].decode() for x in blob.splitlines() if x.startswith(b">")]
        if set(ids)!=set(lengths) or len(ids)!=len(lengths):
            raise ValueError("unexpected FASTA identifiers")
        layouts[name]=ids,hashlib.sha256(blob).hexdigest()
    configs=list(itertools.product(layouts,POLICIES,(1,2)))
    rng=random.Random(91387)
    with (out/"trials.csv").open("x",newline="") as f, (out/"rank_work.csv").open("x",newline="") as g:
        trials=csv.DictWriter(f,TRIAL_FIELDS)
        workers=csv.DictWriter(g,WORK_FIELDS)
        trials.writeheader(); workers.writeheader()
        for block in range(1,6):
            rng.shuffle(configs)
            for order,(layout,policy,nodes) in enumerate(configs,1):
                ids,sha=layouts[layout]
                args=[binary,"matrix","--input",inputs/(layout+".fa"),"--schedule",policy,
                      "--grain","4","--threads","8","--warmups","1","--trials","1","--report-rss","1"]
                result=run_rank_zero(launch(nodes,24,8),args,work)
                check_output(result.stdout,ids,expected)
                info=[json.loads(line) for line in result.stderr.splitlines() if line.startswith("{")]
                rank_data=[row for row in info if "rank" in row]
                [timing]=[row for row in info if "mode" in row]
                ranks=24*nodes; cores=192*nodes
                if {r["rank"] for r in rank_data}!=set(range(ranks)) or len(rank_data)!=ranks:
                    raise ValueError("missing rank diagnostics")
                if (timing["ranks"],timing["threads_per_rank"],timing["sampled_physical_cores"])!=(ranks,8,cores):
                    raise ValueError("unexpected placement")
                if sum(r["pairs"] for r in rank_data)!=len(expected) or sum(r["cells"] for r in rank_data)!=cells:
                    raise ValueError("pair assignment is incomplete")
                src=dict(job=os.environ["SLURM_JOB_ID"],block=block,layout=layout,schedule=policy,nodes=nodes)
                for row in rank_data:
                    workers.writerow({**src,**{k:row[k] for k in WORK_FIELDS if k in row}})
                t=timing["seconds"]
                total_stolen=sum(r["stolen_pairs"] for r in rank_data)
                across=sum(r["cross_node_pairs"] for r in rank_data)
                busy=sum(r["busy_seconds"] for r in rank_data)/(cores*t)
                if not 0<busy<1.001 or not 0<=across<=total_stolen<=len(expected) or (policy!="steal" and total_stolen):
                    raise ValueError("invalid scheduler instrumentation")
                trials.writerow(dict(**src,order=order,ranks=ranks,threads=8,cores=cores,grain=4,pairs=len(expected),
                                     cells=cells,seconds=t,gcups=cells/t/1e9,busy_fraction=busy,
                                     max_mean_rank_cells=max(r["cells"] for r in rank_data)*ranks/cells,
                                     stolen_pairs=total_stolen,cross_node_pairs=across,
                                     max_rank_rss_kib=timing["max_rank_rss_kib"],input_sha256=sha,
                                     canonical_scores_sha256=canonical))
                f.flush(); g.flush()
            print(f"verified block {block} of 5",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("binary","inputs","out","work","windows"):
        parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    measure(*(p.resolve() for p in vars(args).values()))
