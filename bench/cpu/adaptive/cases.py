"""Benchmark inputs with independently checked anchors on the two chromosome sequences."""
import csv
import hashlib
import itertools
import random
from pathlib import Path


def read_pair(path):
    with path.open("rb") as stream:
        return [b"".join(x.strip() for x in group) for header, group in
                itertools.groupby(stream, lambda line: line.startswith(b">")) if not header]


def make_cases(tmp):
    rows = []
    def add(name, path, origin, coords=""):
        a,b=read_pair(path)
        rows.append(dict(case=name,path=path,origin=origin,coordinates=coords,m=len(a),n=len(b),
                         sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    def write(name, a, b, origin, coords=""):
        p=tmp/(name+".fa")
        p.write_bytes(b">A\n"+a+b"\n>B\n"+b+b"\n")
        add(name,p,origin,coords)
    base=Path("data/synthetic")
    with (base/"manifest.tsv").open() as stream:
        for row in csv.DictReader(stream,delimiter="\t"):
            if row["length"]=="20000" and row["rep"]=="1" and float(row["rate"]) in (0.02,0.10):
                add("syn-"+str(round(float(row["rate"])*100)),base/row["path"],"project synthetic")
    rng=random.Random(260128)
    draw=lambda n: bytes(rng.choices(b"ACGT",k=n))
    write("independent",draw(20000),draw(20000),"independent uniform strings, seed 260128")
    full=Path("data/human/chr21-4m.fa")
    assert hashlib.sha256(full.read_bytes()).hexdigest()=="862d699ddbce5424e0583108f31f287213727ed50e1a91b89d27fa649fb35a25"
    a,b=read_pair(full)
    x0,y0=34500000,32882567
    for x1,y1 in ((34525000,32907565),(34600000,32982577),(34750000,33132584),(35500000,33882904)):
        # Coordinates from the positive-strand UCSC chain 23. Check exact unique boundary contexts.
        for x,y in ((x0,y0),(x1,y1)):
            i,j=x-34000000,y-32381795
            flank=a[i-20:i+21]
            assert len(flank)==41 and flank==b[j-20:j+21] and a.count(flank)==b.count(flank)==1
        length=x1-x0
        write(f"chr21-{length//1000}k",a[x0-34000000:x1-34000000],b[y0-32381795:y1-32381795],
              "NC_000021.9 vs NC_060945.1",f"{x0}-{x1}/{y0}-{y1}")
    assert len(rows)==7
    return rows
