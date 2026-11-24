#!/usr/bin/env python3
"""Build two homologous chr21 pairs and a whole-chromosome sparse-edit control."""
import argparse
import csv
import gzip
import hashlib
import shutil
import tempfile
import urllib.request
from pathlib import Path

REFS = {
    "hg38": ("NC_000021.9", "ff47a5bff4a148bc684e444e4c52eaceec1290cb6531b2d8c613becdebbad23c"),
    "hs1": ("NC_060945.1", "f210b99cc28cb8bf66a9aca6d0cd6ceaab1ef8a78ff135876f25644147e27481"),
}
CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
CHAIN_HASH = "926d4583eda8f90c88434293859b06906573f2469630e898b8b475ac8751327a"
START = 34_000_000
ENDS = (36_000_000, 38_000_000)


def download(url, path):
    if path.is_file():
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with urllib.request.urlopen(url, timeout=120) as source:
                shutil.copyfileobj(source, tmp)
            tmp.flush()
            tmp_path.replace(path)
        finally:
            tmp_path.unlink(missing_ok=True)


def chromosome(cache, alias):
    accession, digest = REFS[alias]
    path = cache / (alias + ".fa")
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
           f"&id={accession}&rettype=fasta&retmode=text")
    download(url, path)
    with path.open("rb") as fasta:
        if fasta.readline().split()[0] != (">" + accession).encode():
            raise ValueError(f"unexpected FASTA accession for {alias}")
        seq = b"".join(line.strip() for line in fasta).upper()
    if hashlib.sha256(seq).hexdigest() != digest:
        raise ValueError(f"reference digest mismatch for {accession}")
    return seq


def anchor_points(cache):
    path = cache / "hg38ToHs1.over.chain.gz"
    download(CHAIN, path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != CHAIN_HASH:
        raise ValueError("chain digest mismatch")
    matches = {x: [] for x in (START, *ENDS)}
    active = False
    with gzip.open(path, "rt") as chains:
        for line in chains:
            if line.startswith("chain "):
                _, _, name, _, strand, start, end, query, _, qstrand, qstart, _, cid = line.split()
                t, q = int(start), int(qstart)
                active = name == "chr21" and strand == "+" and any(t <= x < int(end) for x in matches)
            elif active and line.strip():
                fields = list(map(int, line.split()))
                size = fields[0]
                for x in matches:
                    if t <= x < t + size:
                        matches[x].append((cid, query, qstrand, q + x - t))
                if len(fields) == 3:
                    t += size + fields[1]
                    q += size + fields[2]
    for x, hits in matches.items():
        if len(hits) != 1 or hits[0][:3] != ("23", "chr21", "+"):
            raise ValueError(f"ambiguous or unexpected chain at hg38 chr21:{x}")
    return {x: hits[0][3] for x, hits in matches.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("data/human"))
    parser.add_argument("--manifest", type=Path, default=Path("results/human/manifest.tsv"))
    parser.add_argument("--cache", type=Path, help="reuse downloaded reference files")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        cache = args.cache or Path(tmp)
        cache.mkdir(parents=True, exist_ok=True)
        hg, hs = (chromosome(cache, alias) for alias in ("hg38", "hs1"))
        points = anchor_points(cache)

    rows = []

    def pair(name, accession_a, start_a, end_a, a, accession_b, start_b, end_b, b, edit):
        path = args.out / (name + ".fa")
        with path.open("wb") as stream:
            for tag, sequence in (("A", a), ("B", b)):
                stream.write(f">{name}_{tag}\n".encode())
                for i in range(0, len(sequence), 80):
                    stream.write(sequence[i:i+80] + b"\n")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(dict(case=name, accession_a=accession_a, start_a=start_a, end_a=end_a,
                         accession_b=accession_b, start_b=start_b, end_b=end_b, edit=edit,
                         m=len(a), n=len(b), sha256=digest))

    hg_id, hs_id = REFS["hg38"][0], REFS["hs1"][0]
    for end in ENDS:
        a, b = hg[START:end], hs[points[START]:points[end]]
        for seq, window in ((hg, a), (hs, b)):
            if seq.count(window[:41]) != 1 or seq.count(window[-41:]) != 1 or b"N" in window:
                raise ValueError("window flank is nonunique or sequence contains N")
        pair(f"chr21-{(end - START) // 1_000_000}m", hg_id, START, end, a,
             hs_id, points[START], points[end], b, "public assemblies, chain 23 (+)")
    cuts = [0, *(k * len(hs) // 65 for k in range(1, 65)), len(hs)]
    mutant = b"G".join(hs[left:right] for left, right in zip(cuts, cuts[1:]))
    pair("chr21-full-insert64", hs_id, 0, len(hs), hs, hs_id, 0, len(hs), mutant,
         "insert G at 64 evenly spaced reference boundaries")
    with args.manifest.open("w", newline="") as output:
        writer = csv.DictWriter(output, rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
