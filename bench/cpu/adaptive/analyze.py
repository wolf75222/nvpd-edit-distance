#!/usr/bin/env python3
"""Rebuild the statistics of the adaptive CPU study from per-process trials and counters."""
import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w",newline="") as stream:
        out=csv.DictWriter(stream,rows[0].keys());out.writeheader();out.writerows(rows)


def summary(values):
    assert len(values)==5 and all(math.isfinite(v) for v in values)
    mean,sd=statistics.mean(values),statistics.stdev(values)
    half=2.7764451052*sd/math.sqrt(5)
    return dict(n=5,mean=mean,sd=sd,low95=mean-half,high95=mean+half)


def analyze(folder):
    inputs={r["case"]:r for r in read(folder/"cases.csv")}
    trials=read(folder/"trials.csv")
    timed=defaultdict(list)
    for r in trials:
        assert r["distance"]==inputs[r["case"]]["distance"]
        timed[r["case"],r["algorithm"]].append(r)
    rows=[];work=[]
    for (name,alg),g in sorted(timed.items()):
        assert len({v["block"] for v in g})==5
        times=[float(v["seconds"]) for v in g]
        rows.append(dict(case=name,algorithm=alg,metric="seconds",**summary(times)))
        if alg!="nvpd-adaptive":continue
        assert len({(v["cells"],v["threshold"],v["passes"],v["prefix"],v["suffix"]) for v in g})==1
        r=g[0];case=inputs[name]
        m,n=sorted([int(case["m"]),int(case["n"])])
        clip=int(r["prefix"])+int(r["suffix"]);k=int(r["threshold"])
        delta=n-m;width=delta+2*((k-delta)//2)+1
        cells=int(r["cells"])
        work.append(dict(case=name,m=m,n=n,distance=r["distance"],cells=cells,
                         passes=r["passes"],threshold=k,prefix=r["prefix"],suffix=r["suffix"],
                         cell_ratio=m*n/cells,visited_fraction=cells/(m*n),
                         dense_array_bytes=4*6*(n+1),
                         adaptive_array_bytes=4*(4*(n-clip+1)+2*(width+2))))
        for baseline in ("nvpd-omp","band-control","banded-trim"):
            reference=timed.get((name,baseline),[])
            if not reference:continue
            paired={v["block"]:float(v["seconds"]) for v in reference}
            ratio=[paired[v["block"]]/float(v["seconds"]) for v in g]
            rows.append(dict(case=name,algorithm=baseline,metric="paired_time_over_adaptive",
                             **summary(ratio)))
    write(folder/"summary.csv",rows);write(folder/"work.csv",work)

    grouped=defaultdict(list)
    checksums=defaultdict(set)
    for r in read(folder/"calibration.csv"):
        cells=int(r["elements"])*int(r["rounds"]);t=float(r["seconds"])
        kind=r["kernel"];factor={"read":4,"triad":12,"row":1}[kind]
        grouped[kind,int(r["footprint_kib"])].append(factor*cells/t*1e-9)
        checksums[kind,int(r["footprint_kib"])].add((int(r["elements"]), int(r["rounds"]),r["checksum"]))
    assert all(len(values)==1 for values in checksums.values())
    write(folder/"calibration-summary.csv",[dict(kernel=kind,footprint_kib=size,
          unit="Gcell/s" if kind=="row" else "GB/s",**summary(values))
          for (kind,size),values in sorted(grouped.items())])

    pmu=defaultdict(dict);samplemeta={}
    for r in read(folder/"counters.csv"):
        key=(r["case"],r["algorithm"],r["block"]);event=r["event"].removesuffix(":u")
        assert event not in pmu[key] and float(r["running_pct"])>=95
        pmu[key][event]=int(r["count"])
        if event=="r4844":samplemeta[key]=(int(r["cells"]),float(r["seconds"]))
    summary_pmu=defaultdict(list)
    for key,c in pmu.items():
        name,alg,_=key;cells,seconds=samplemeta[key]
        read_bytes=64*(c["r4844"]+c["rff72"])
        assert read_bytes>0
        values=dict(estimated_dram_bytes_per_update=read_bytes/cells,
                    updates_per_estimated_dram_byte=cells/read_bytes,
                    gcell_s=cells/seconds*1e-9)
        if name=="calibration-dram":
            values["estimated_over_known_read_bytes"]=read_bytes/(4*cells)
        else:
            values.update(ipc=c["instructions"]/c["cycles"],cycles_per_cell=c["cycles"]/cells,
                          instructions_per_cell=c["instructions"]/cells,
                          l1_misses_per_1000_cells=1000*c["L1-dcache-load-misses"]/cells,
                          branch_misses_per_1000_cells=1000*c["branch-misses"]/cells)
        for metric,value in values.items():summary_pmu[name,alg,metric].append(value)
    write(folder/"pmu-summary.csv",[dict(case=name,algorithm=alg,metric=metric,**summary(values))
          for (name,alg,metric),values in sorted(summary_pmu.items())])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("result_folder",type=Path)
    analyze(p.parse_args().result_folder)
