"""Compare three exact methods on chromosome windows in five matched blocks."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import os
from pathlib import Path
import re
import resource
import statistics
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2] / "human/algorithms"))
from study import read_pair, replay

CASES = (("chr21-2m",5599,"2b6c4c29b210d67fc56db9e32af26b0610a742804349e6728d94d123b6af877d"),
         ("chr21-4m",12482,"862d699ddbce5424e0583108f31f287213727ed50e1a91b89d27fa649fb35a25"),
         ("chr21q-8m",47493,"330dd197a267a23c63aa9a5f60aa849665923ddaaa6cee5ca6cec9e69627a022"))
METHODS = ("astarpa2-full","quicked-certified","biwfa")
COLUMNS = "case algorithm block order cpu status seconds peak_rss_kib distance source retries cigar_bytes cigar_sha256 input_sha256 reason".split()


def lane(arguments):
    binary, inputs, out, spec, cpu = arguments
    case, expected, fingerprint = spec
    os.sched_setaffinity(0,{cpu})
    if os.sched_getaffinity(0) != {cpu}:
        raise RuntimeError("invalid CPU affinity")
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    path = inputs/(case+".fa")
    if hashlib.sha256(path.read_bytes()).hexdigest() != fingerprint:
        raise ValueError("unexpected FASTA digest")
    a,b = read_pair(path)
    capped = {}
    with (out/(case+".csv")).open("x",newline="") as f:
        writer=csv.DictWriter(f,COLUMNS)
        writer.writeheader()
        for block in range(5):
            order = list(METHODS if block < 3 else METHODS[::-1])
            shift = block % len(order)
            for pos,method in enumerate(order[shift:]+order[:shift],1):
                row=dict(case=case,algorithm=method,block=block+1,order=pos,cpu=cpu,input_sha256=fingerprint,status="ok")
                if method in capped:
                    row.update(status="not_attempted",reason=capped[method])
                else:
                    try:
                        result=subprocess.run([binary,path,method],timeout=120,capture_output=True)
                        if result.returncode:
                            detail = result.stderr.decode(errors="replace").lower()
                            code = "memory_limit" if any(word in detail for word in ("alloc","memory")) else "backend_error"
                            row.update(status=code,reason=f"exit code {result.returncode}")
                            capped[method] = code
                        else:
                            head, cigar = result.stdout.split(b"\n",1)
                            elapsed,rss,score,source,retries = head.decode().split()
                            cigar = cigar.rstrip(b"\n")
                            if int(score)!=expected or float(elapsed)<=0 or int(rss)<=0:
                                raise ValueError(f"invalid result for {case}/{method}")
                            if method == "astarpa2-full":
                                normalized=re.sub(rb"(?<![0-9])([=XID])",rb"1\1",cigar)
                                replay(a,b,normalized,expected)
                            elif method == "biwfa":
                                replay(a,b,cigar,expected)
                            elif source == "edlib-fallback":
                                row["status"] = "score_fallback"
                            elif source not in ("quicked-default","quicked-widened"):
                                raise ValueError(f"unknown QuickEd result source {source}")
                            row.update(seconds=elapsed,peak_rss_kib=rss,distance=score,source=source,retries=retries,
                                       cigar_bytes=len(cigar),cigar_sha256=hashlib.sha256(cigar).hexdigest() if cigar else "")
                    except subprocess.TimeoutExpired:
                        row.update(status="timeout",reason="120 s process limit, includes warmup and measured call")
                        capped[method] = "earlier timeout"
                writer.writerow(row)
                f.flush()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("binary","inputs","out"): p.add_argument(key,type=Path)
    a=p.parse_args()
    cpus=sorted(os.sched_getaffinity(0))
    if len(cpus)<len(CASES): raise ValueError("need three distinct physical cores")
    with ProcessPoolExecutor(max_workers=len(CASES)) as pool:
        list(pool.map(lane,[(a.binary,a.inputs,a.out,case,cpu) for case,cpu in zip(CASES,cpus)]))
    rows=[]
    for name,_,_ in CASES:
        with (a.out/(name+".csv")).open() as f: rows.extend(csv.DictReader(f))
    with (a.out/"trials.csv").open("x",newline="") as f:
        w=csv.DictWriter(f,COLUMNS);w.writeheader();w.writerows(rows)
    with (a.out/"summary.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("case","algorithm","metric","scheduled","completed","mean","sd","low95","high95"))
        for name,_,_ in CASES:
            for method in METHODS:
                valid=[r for r in rows if (r["case"],r["algorithm"],r["status"])==(name,method,"ok")]
                for field in ("seconds","peak_rss_kib"):
                    x=[float(r[field]) for r in valid]
                    if len(x)<5:
                        w.writerow((name,method,field,5,len(x),"","","",""))
                        continue
                    avg,sd=statistics.mean(x),statistics.stdev(x)
                    half=2.7764451051977987*sd/len(x)**0.5
                    w.writerow((name,method,field,5,len(x),avg,sd,avg-half,avg+half))


if __name__=="__main__":
    main()
