"""Extract four chromosome pairs from the two published HG002 BGZF assemblies."""
import argparse
import bisect
import concurrent.futures
import csv
import gzip
import hashlib
import io
from pathlib import Path
import struct
import time
import urllib.request

BASE = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies"
VCF = ("https://raw.githubusercontent.com/marbl/HG002-issues/44beefe41da04891c5e68b16fc1545ca82522693/"
       "patches/already_applied/hg002v1.1_to_hg002v1.2.vcf.gz")
CHROMS = [f"chr{chrom}_{hap}" for chrom in (1, 21) for hap in ("MATERNAL", "PATERNAL")]


def fetch(url, interval=None):
    headers = {} if interval is None else {"Range": "bytes=%d-%d" % interval}
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as reply:
                if interval and (reply.status != 206 or
                                 not reply.headers["Content-Range"].startswith("bytes %d-%d/" % interval)):
                    raise ValueError(f"unexpected byte range: {reply.headers.get('Content-Range')}")
                return reply.read()
        except OSError:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def extract(version, out):
    url = f"{BASE}/hg002v{version}.fasta.gz"
    fai = {r[0]: tuple(map(int, r[1:5])) for r in
           (line.split() for line in fetch(url + ".fai").decode().splitlines())}
    index = fetch(url + ".gzi")
    offsets = [(0, 0)] + list(struct.iter_unpack("<QQ", index[8:]))
    assert len(offsets) - 1 == struct.unpack("<Q", index[:8])[0]
    uncompressed = [u for _, u in offsets]
    for chrom in CHROMS:
        length, offset, bases, width = fai[chrom]
        end = offset + ((length - 1) // bases) * width + (length - 1) % bases + 1
        lo = bisect.bisect_right(uncompressed, offset) - 1
        hi = bisect.bisect_left(uncompressed, end)
        a, b = offsets[lo][0], offsets[hi][0]
        chunks = [(pos, min(pos + 4 * 1024**2, b) - 1) for pos in range(a, b, 4 * 1024**2)]
        with concurrent.futures.ThreadPoolExecutor(6) as pool:
            compressed = b"".join(pool.map(lambda interval: fetch(url, interval), chunks))
        start = offset - offsets[lo][1]
        stream = gzip.GzipFile(fileobj=io.BytesIO(compressed))
        stream.seek(start)
        record = stream.read(end - offset)
        seq = b"".join(record.splitlines())
        if len(seq) != length:
            raise ValueError(f"{chrom} length differs from the FASTA index")
        path = out / f"{chrom}.v{version}.fa"
        with path.open("wb") as f:
            f.write(f">HG002v{version}_{chrom}\n".encode())
            f.write(record)
            f.write(b"\n")
        print(f"{chrom} {version} {length}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    for version in ("1.1", "1.2"):
        extract(version, out)
    patches = {chrom: [] for chrom in CHROMS}
    for line in gzip.decompress(fetch(VCF)).splitlines():
        if line.startswith(b"#"):
            continue
        chrom, pos, _, ref, alt = line.split(b"\t")[:5]
        if chrom.decode() in patches:
            patches[chrom.decode()].append((int(pos) - 1, ref, alt))
    with (out / "manifest.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(("chrom", "version", "bases", "sha256_bases", "n_ambiguous", "vcf_records"))
        for chrom in CHROMS:
            seqs = []
            for version in ("1.1", "1.2"):
                with (out / f"{chrom}.v{version}.fa").open("rb") as fasta:
                    next(fasta)
                    seq = fasta.read().replace(b"\n", b"")
                if seq.translate(None, b"ACGTN"):
                    raise ValueError(f"unexpected base in {chrom}")
                writer.writerow((chrom, version, len(seq), hashlib.sha256(seq).hexdigest(),
                                 seq.count(b"N"), len(patches[chrom])))
                seqs.append(seq)
            pos = 0
            digest = hashlib.sha256()
            for i, ref, alt in patches[chrom]:
                if i < pos or seqs[0][i:i + len(ref)] != ref:
                    raise ValueError(f"inconsistent VCF record in {chrom} at {i + 1}")
                digest.update(memoryview(seqs[0])[pos:i])
                digest.update(alt)
                pos = i + len(ref)
            digest.update(memoryview(seqs[0])[pos:])
            if digest.digest() != hashlib.sha256(seqs[1]).digest():
                raise ValueError(f"patched v1.1 differs from v1.2 in {chrom}")


if __name__ == "__main__":
    main()
