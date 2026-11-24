"""Aggregate independent launches and within-block one-to-two-node ratios."""
import argparse
import bisect
import csv
import itertools
import math
from pathlib import Path
import statistics

METRICS = ("seconds", "gcups", "busy_fraction", "max_mean_rank_cells",
           "stolen_pairs", "cross_node_pairs", "max_rank_rss_kib")
T4 = 2.7764451051977987


def estimate(values):
    if len(values)!=5 or not all(math.isfinite(x) for x in values):
        raise ValueError("expected five complete independent blocks")
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    half = T4 * sd / math.sqrt(len(values))
    return dict(n=len(values),mean=mean,sd=sd,low95=mean-half,high95=mean+half)


def save(path, records):
    with path.open("w",newline="") as stream:
        writer = csv.DictWriter(stream,records[0].keys())
        writer.writeheader()
        writer.writerows(records)


def predicted_loads(lengths,policy,ranks):
    """Invert the Hilbert index independently, then predict rank work in exact integers."""
    n=len(lengths)
    pairs=list(itertools.combinations(range(n),2))
    if policy.startswith("hilbert"):
        side=1 << (n-1).bit_length()
        pairs=[]
        for code in range(side*side):
            x=y=0
            digit=code
            scale=1
            while scale<side:
                rx=(digit//2)&1
                ry=(digit^rx)&1
                if ry==0:
                    if rx:
                        x,y=scale-1-x,scale-1-y
                    x,y=y,x
                x+=scale*rx
                y+=scale*ry
                digit//=4
                scale*=2
            if x<y<n:
                pairs.append((x,y))
    weights=[lengths[i]*lengths[j] for i,j in pairs]
    if policy=="cyclic":
        return [(len(weights[r::ranks]),sum(weights[r::ranks])) for r in range(ranks)]
    prefix=[0,*itertools.accumulate(weights)]
    cuts=[i*len(weights)//ranks for i in range(ranks+1)]
    if policy=="hilbert-weighted" and prefix[-1]:
        for r in range(1,ranks):
            target=prefix[-1]*r//ranks
            hi=bisect.bisect_left(prefix,target)
            cuts[r]=hi-1 if hi and target-prefix[hi-1]<prefix[hi]-target else hi
    return [(b-a,prefix[b]-prefix[a]) for a,b in zip(cuts,cuts[1:])]


def summarize(path):
    with (path/"trials.csv").open() as f:
        rows=list(csv.DictReader(f))
    with (path/"rank_work.csv").open() as f:
        work=list(csv.DictReader(f))
    choices = ("static","cyclic","hilbert","steal","hilbert-weighted")
    policies = [s for s in choices if s in {r["schedule"] for r in rows}]
    if len(policies) not in (4,5) or len(rows)!=20*len(policies) or len(work)!=720*len(policies):
        raise ValueError("incomplete human window study")
    grouped={}
    for r in rows:
        tag=r["layout"],r["schedule"],int(r["nodes"]),int(r["block"])
        if tag in grouped: raise ValueError("duplicate launch")
        grouped[tag]=r
    detail={}
    for r in work:
        tag=r["layout"],r["schedule"],int(r["nodes"]),int(r["block"])
        detail.setdefault(tag,[]).append(r)
    with (path.parent/"windows.tsv").open() as f:
        windows=list(csv.DictReader(f,delimiter="\t"))
    order={"locus":windows,"length":sorted(windows,key=lambda w:(int(w["length"]),w["id"]))}
    loads={(layout,policy,nodes):predicted_loads([int(x["length"]) for x in order[layout]],policy,24*nodes)
           for layout,policy,nodes in itertools.product(order,policies,(1,2)) if policy!="steal"}
    for key,r in grouped.items():
        local=detail[key]
        if {int(x["rank"]) for x in local}!=set(range(int(r["ranks"]))) or len(local)!=int(r["ranks"]):
            raise ValueError(f"invalid ranks: {key}")
        for total,field in (("pairs","pairs"),("cells","cells"),("stolen_pairs","stolen_pairs"),
                            ("cross_node_pairs","cross_node_pairs")):
            if int(r[total])!=sum(int(x[field]) for x in local):
                raise ValueError(f"inconsistent {field}: {key}")
        if r["schedule"]!="steal":
            observed=[(int(x["pairs"]),int(x["cells"])) for x in sorted(local,key=lambda x:int(x["rank"]))]
            if observed!=loads[key[:3]]:
                raise ValueError(f"rank counts or work differ from the independent partition: {key}")
        expected = sum(float(x["busy_seconds"]) for x in local)/(float(r["seconds"])*int(r["cores"]))
        if not math.isclose(expected,float(r["busy_fraction"]),rel_tol=1e-10):
            raise ValueError("inconsistent busy fraction")
    summaries=[]
    scaling=[]
    for layout,policy in itertools.product(("locus","length"),policies):
        conditions=dict(layout=layout,schedule=policy)
        for nodes in (1,2):
            series=[grouped[layout,policy,nodes,b] for b in range(1,6)]
            for metric in METRICS:
                summaries.append(dict(**conditions,nodes=nodes,metric=metric,
                                      **estimate([float(r[metric]) for r in series])))
        gain=[float(grouped[layout,policy,1,b]["seconds"])/float(grouped[layout,policy,2,b]["seconds"])
              for b in range(1,6)]
        for metric,values in (("speedup",gain),("efficiency",[x/2 for x in gain])):
            scaling.append(dict(**conditions,metric=metric,**estimate(values)))
    comparisons = []
    contenders = [("static","cyclic"), ("static","steal"), ("static","hilbert")]
    if "hilbert-weighted" in policies:
        contenders += [("hilbert","hilbert-weighted"),("cyclic","hilbert-weighted"),
                       ("steal","hilbert-weighted")]
    for layout,nodes,(baseline,candidate) in itertools.product(("locus","length"),(1,2),contenders):
        ratios = [float(grouped[layout,baseline,nodes,b]["seconds"])
                  /float(grouped[layout,candidate,nodes,b]["seconds"]) for b in range(1,6)]
        comparisons.append(dict(layout=layout,nodes=nodes,baseline=baseline,candidate=candidate,
                                metric="time_ratio",**estimate(ratios)))
    save(path/"paired.csv",comparisons)
    save(path/"summary.csv",summaries)
    save(path/"scaling.csv",scaling)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path",type=Path)
    summarize(parser.parse_args().path)
