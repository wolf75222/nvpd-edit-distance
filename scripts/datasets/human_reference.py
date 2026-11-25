#!/usr/bin/env python3
"""Prepare homologous chromosome windows and whole human chromosome pairs."""
import argparse
import csv
import gzip
import hashlib
import io
import shutil
import tempfile
import urllib.request
from pathlib import Path

REFS = {
    "chm20v1": ("CP068258.1", "5b9b1e5d051827ad0b1e693d6868bf43edcbba553187aa09c1034decb2eaeb1e"),
    "chm20v2": ("CP068258.2", "a115f944c3e4753444d189a5518e32c086bbd73e052a3f227f2a6ae6211f584f"),
    "hg38": ("NC_000021.9", "ff47a5bff4a148bc684e444e4c52eaceec1290cb6531b2d8c613becdebbad23c"),
    "hs1": ("NC_060945.1", "f210b99cc28cb8bf66a9aca6d0cd6ceaab1ef8a78ff135876f25644147e27481"),
}
CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
CHAIN_HASH = "926d4583eda8f90c88434293859b06906573f2469630e898b8b475ac8751327a"
START = 34_000_000
ENDS = (36_000_000, 38_000_000)


def verified_cache(url, path, read):
    """Publish a download after its reader has verified the pinned content."""
    if path.is_file():
        try:
            return read(path)
        except ValueError:
            pass  # Re-fetch a damaged cached copy.
    output = tempfile.NamedTemporaryFile(dir=path.parent, delete=False)
    tmp_path = Path(output.name)
    try:
        with output, urllib.request.urlopen(url, timeout=120) as source:
            shutil.copyfileobj(source, output)
        value = read(tmp_path)
        tmp_path.replace(path)
        return value
    finally:
        tmp_path.unlink(missing_ok=True)


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def chromosome(cache, alias):
    accession, digest = REFS[alias]
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
           f"&id={accession}&rettype=fasta&retmode=text")

    def read(path):
        with path.open("rb") as fasta, io.BytesIO() as bases:
            if fasta.readline().split()[:1] != [(">" + accession).encode("ascii")]:
                raise ValueError(f"unexpected FASTA accession for {alias}")
            for line in fasta:
                bases.write(line.strip().upper())
            seq = bases.getvalue()
        if hashlib.sha256(seq).hexdigest() != digest:
            raise ValueError(f"reference digest mismatch for {accession}")
        return seq

    return verified_cache(url, cache / (alias + ".fa"), read)


def chain_path(cache):
    path = cache / "hg38ToHs1.over.chain.gz"

    def verify(candidate):
        if file_sha256(candidate) != CHAIN_HASH:
            raise ValueError("UCSC chain digest mismatch")

    verified_cache(CHAIN, path, verify)
    return path


def alignment_blocks(cache):
    """Aligned intervals (hg38 start, hg38 end, hs1 start) on chr21 chain 23."""
    path = chain_path(cache)
    blocks = []
    active = False
    with gzip.open(path, "rt") as source:
        for line in source:
            if line.startswith("chain "):
                _, _, chrom, _, strand, pos, _, query, _, sense, qpos, _, chain = line.split()
                active = (chrom, strand, query, sense, chain) == ("chr21", "+", "chr21", "+", "23")
                x, y = int(pos), int(qpos)
            elif active and line.strip():
                fields = [int(v) for v in line.split()]
                blocks.append((x, x + fields[0], y))
                if len(fields) == 3:
                    x += fields[0] + fields[1]
                    y += fields[0] + fields[2]
    if not blocks:
        raise ValueError("chain 23 is absent")
    return blocks


def anchor_points(cache):
    path = chain_path(cache)
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
    parser.add_argument("--whole-assemblies", action="store_true", help="also compare complete chromosome assemblies")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        cache = args.cache or Path(tmp)
        cache.mkdir(parents=True, exist_ok=True)
        hg, hs = (chromosome(cache, alias) for alias in ("hg38", "hs1"))
        points = anchor_points(cache)
        if args.whole_assemblies:
            v1, v2 = (chromosome(cache, key) for key in ("chm20v1", "chm20v2"))

    rows = []

    def pair(name, accession_a, start_a, end_a, a, accession_b, start_b, end_b, b, edit):
        path = args.out / (name + ".fa")
        with path.open("wb") as stream:
            for tag, sequence in (("A", a), ("B", b)):
                stream.write(f">{name}_{tag}\n".encode())
                for i in range(0, len(sequence), 80):
                    stream.write(sequence[i:i+80] + b"\n")
        digest = file_sha256(path)
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
    if args.whole_assemblies:
        pair("chr20-assemblies", REFS["chm20v1"][0], 0, len(v1), v1, REFS["chm20v2"][0],
             0, len(v2), v2, "T2T-CHM13 assemblies v1.0 and v2.0")
        pair("chr21-assemblies", hg_id, 0, len(hg), hg, hs_id, 0, len(hs), hs,
             "complete GRCh38.p14 and T2T-CHM13v2.0")
    with args.manifest.open("w", newline="") as output:
        writer = csv.DictWriter(output, rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
