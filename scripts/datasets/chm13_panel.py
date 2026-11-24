#!/usr/bin/env python3
"""Prepare the non-acrocentric chromosomes of two public T2T-CHM13 assemblies."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen

SOURCE = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies"
RELEASES = {"1.0": "6d827b6512562630137008830c46e1ac", "1.1": "1cab2b2776005cdf339ec9f283ba2c70"}
CHROMOSOMES = tuple(f"chr{i}" for i in range(1, 23) if i not in (13, 14, 15, 21, 22)) + ("chrX",)


def digest(path, algorithm):
    result = hashlib.new(algorithm)
    with path.open("rb") as reader:
        for block in iter(lambda: reader.read(4 << 20), b""):
            result.update(block)
    return result.hexdigest()


def archive(cache, version):
    path = cache / f"chm13.v{version}.fa.gz"
    if not path.exists():
        url = f"{SOURCE}/chm13.draft_v{version}.fasta.gz"
        with tempfile.NamedTemporaryFile(dir=cache, delete=False) as stream:
            tmp = Path(stream.name)
        try:
            with urlopen(url, timeout=120) as source, tmp.open("wb") as sink:
                shutil.copyfileobj(source, sink, 4 << 20)
            if digest(tmp, "md5") != RELEASES[version]:
                raise ValueError(f"archive checksum mismatch: {url}")
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
    elif digest(path, "md5") != RELEASES[version]:
        raise ValueError(f"cached archive checksum mismatch: {path}")
    return path


def extract(archive_path, version, folder):
    lengths = {}
    output = None
    try:
        with gzip.open(archive_path, "rb") as genome:
            for line in genome:
                if line.startswith(b">"):
                    if output is not None:
                        output.close()
                        output = None
                    name = line[1:].split()[0].decode("ascii")
                    if name in CHROMOSOMES:
                        if name in lengths:
                            raise ValueError(f"duplicate {name} in release {version}")
                        lengths[name] = 0
                        output = (folder / f"{name}.fa").open("ab")
                        output.write(f">CHM13v{version}_{name}\n".encode("ascii"))
                elif output is not None:
                    bases = line.rstrip().upper()
                    lengths[name] += len(bases)
                    output.write(bases + b"\n")
    finally:
        if output is not None:
            output.close()
    if set(lengths) != set(CHROMOSOMES):
        raise ValueError(f"incomplete panel in release {version}")
    return lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("data/human/chm13-panel"))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out) as temporary:
        folder = Path(temporary)
        sizes = [extract(archive(args.cache, version), version, folder) for version in RELEASES]
        with (folder / "manifest.tsv").open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(("case", "m", "n", "bytes", "input_sha256"))
            for name in CHROMOSOMES:
                path = folder / f"{name}.fa"
                writer.writerow((name, sizes[0][name], sizes[1][name], path.stat().st_size,
                                 digest(path, "sha256")))
        for name in CHROMOSOMES:
            (folder / f"{name}.fa").replace(args.out / f"{name}.fa")
        (folder / "manifest.tsv").replace(args.out / "manifest.tsv")


if __name__ == "__main__":
    main()
