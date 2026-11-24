import csv
import os
from pathlib import Path
import random
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import measure as protocol


def main():
    data, build, out = map(Path, sys.argv[1:])
    lane = int(os.environ["SLURM_PROCID"])
    name = protocol.READ_CASES[lane]
    case = next(r for r in csv.DictReader((data / "manifest.tsv").open(), delimiter="\t") if r["case"]==name)
    path = data / (name + ".fa")
    cores = sorted(os.sched_getaffinity(0))
    if len(cores) != 1: raise RuntimeError(cores)
    common = dict(job=os.environ["SLURM_JOB_ID"], case=name, cpu=cores[0])
    rng = random.Random(710 + lane)
    variants = ["portable", "blocks"]
    # Verify both builds against an independently decoded complete chromosome pair.
    expected = protocol.canonical_file_digest(path)
    for variant in variants:
        if protocol.normalized_digest(build / variant, path) != expected:
            raise ValueError((name, variant))
    with (out/f"blocks-{lane}.csv").open("x") as sink:
        report = protocol.writer(sink, "job case cpu method block warm_seconds seconds peak_rss_kib instructions cycles L1_loads L1_misses allocations requested_bytes")
        for block in range(1,6):
            rng.shuffle(variants)
            for variant in variants:
                binary=build/variant
                time=protocol.run_trial([binary,path,"reader"],case)
                del time["distance"]
                lines=subprocess.check_output([binary,path,"profile"],text=True,timeout=180).splitlines()[2:]
                values={}
                for line in lines:
                    event,value,*pct=line.split(",")
                    if pct and float(pct[0]) < 99: raise ValueError((event,pct))
                    values[event]=int(value)
                values["L1_loads"]=values.pop("L1-dcache-loads")
                values["L1_misses"]=values.pop("L1-dcache-load-misses")
                report.writerow(dict(common, method=variant, block=block, **time, **values))
                sink.flush()


if __name__=="__main__":
    main()
