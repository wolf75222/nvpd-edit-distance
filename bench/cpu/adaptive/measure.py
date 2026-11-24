#!/usr/bin/env python3
import argparse
import csv
import itertools
import json
import os
import random
import resource
import subprocess
from pathlib import Path
from cases import make_cases


def write(path, rows):
    with path.open("w", newline="") as stream:
        out=csv.DictWriter(stream,rows[0].keys())
        out.writeheader()
        out.writerows(rows)


def cap():
    resource.setrlimit(resource.RLIMIT_AS,(6<<30,6<<30))


def run(cmd, env=None):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, env=env,
                          preexec_fn=cap, timeout=180)


def main():
    p=argparse.ArgumentParser()
    for key in ("build","tmp","out","perf"):
        p.add_argument(key,type=Path)
    args=p.parse_args()
    cli=args.build/"editdist"
    driver=args.tmp/"driver"
    cases=make_cases(args.tmp)
    meta={r["case"]:r for r in cases}
    for c in cases:
        cmd=[cli,"pair","--algorithm","edlib","--input",c["path"]]
        c["distance"]=json.loads(run(cmd).stdout)["distance"]
    write(args.out/"cases.csv",[{k:v for k,v in row.items() if k!="path"} for row in cases])
    rng=random.Random(260128)
    variants=("nvpd-adaptive","band-control","nvpd-omp","banded-trim","edlib","wfa2","biwfa")
    configs=[(name,alg) for name,alg in itertools.product(meta,variants)
             if not (name=="chr21-1000k" and alg=="nvpd-omp")]
    def align(argv, name, alg, env=None):
        c=meta[name]
        exe = args.tmp/"driver-control" if alg=="band-control" else driver
        method = "nvpd-adaptive" if alg=="band-control" else alg
        output=run([*argv,exe,"align",c["path"],method],env)
        t,d,cells,passes,k,prefix,suffix=output.stdout.strip().split(",")
        assert int(d)==c["distance"],(name,alg,d,c["distance"])
        return dict(case=name,algorithm=alg,seconds=float(t),distance=d,cells=cells,
                    passes=passes,threshold=k,prefix=prefix,suffix=suffix)
    # Count events only in the second backend call. Separate passes fit the Zen 4 PMU.
    ctl,ack=args.tmp/"ctl",args.tmp/"ack"
    os.mkfifo(ctl);os.mkfifo(ack)
    env=dict(os.environ,CHPS_PERF_CTL=str(ctl),CHPS_PERF_ACK=str(ack),LC_ALL="C")
    PMU=(("cycles:u","instructions:u","branches:u","branch-misses:u"),
         ("L1-dcache-loads:u","L1-dcache-load-misses:u","r4844:u","rff72:u"))
    def counter_cmd(events):
        return [args.perf,"stat","-x",";","--no-big-num","--delay=-1","--control",
                f"fifo:{ctl},{ack}","-o",args.tmp/"perf.csv","-e",",".join(events),"--"]
    def values():
        for line in (args.tmp/"perf.csv").read_text().splitlines():
            fields=[item.strip() for item in line.split(";")]
            if len(fields)<5 or not fields[2]:continue
            if not fields[0].isdigit() or float(fields[4])<95:
                raise ValueError(f"unsupported or multiplexed event: {line}")
            yield dict(event=fields[2],count=int(fields[0]),running_pct=float(fields[4]))
    counters=[];calib=[];timings=[]
    for block in range(1,6):
        series=[(x,kib) for x in ("read","triad") for kib in (12,384,12288,196608)]+[("row",12)]
        rng.shuffle(series)
        for kind,kib in series:
            cmd=[driver,"calibrate",kind,str(kib)]
            t,n,r,checksum=run(cmd).stdout.strip().split(",")
            assert int(checksum)>0
            calib.append(dict(block=block,kernel=kind,footprint_kib=kib,seconds=float(t),
                              elements=int(n),rounds=int(r),checksum=checksum))
            if kind=="read" and kib==196608:
                args2=[*counter_cmd(PMU[1]),*cmd]
                elapsed,_,_,check=run(args2,env).stdout.strip().split(",")
                assert checksum==check
                for v in values():
                    counters.append(dict(block=block,case="calibration-dram",algorithm="read",seconds=elapsed,
                                         cells=int(n)*int(r),**v))
        write(args.out/"calibration.csv",calib)
        for name in ("syn-10","chr21-100k","independent"):
            for alg in ("nvpd-adaptive","band-control","nvpd-omp"):
                for events in PMU:
                    result=align(counter_cmd(events),name,alg,env)
                    for v in values():
                        counters.append(dict(block=block,case=name,algorithm=alg,seconds=result["seconds"],
                                             cells=result["cells"],**v))
        write(args.out/"counters.csv",counters)
        order=configs.copy()
        rng.shuffle(order)
        for pos,(name,alg) in enumerate(order,1):
            timings.append(dict(block=block,order=pos,**align([],name,alg)))
            write(args.out/"trials.csv",timings)


if __name__=="__main__":
    main()
