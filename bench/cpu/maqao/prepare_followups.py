"""Materialize four changes to the blocked SIMD kernel and the optional layout candidates."""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpu", type=Path, help="src/cpu from the checkout of 50d9bd9")
    parser.add_argument("destination", type=Path, help="new directory for candidate source trees")
    parser.add_argument("--layout", action="store_true", help="also prepare the 40,000-base layout candidates")
    args = parser.parse_args()
    header = args.cpu / "nvpd_avx512.hpp"
    digest = "007caa314c131a8b47e0daaaf195a6d589fc0b9e86a1ec87d190e1a8529c723b"
    if hashlib.sha256(header.read_bytes()).hexdigest() != digest:
        raise ValueError("unexpected reference kernel; use the checkout of 50d9bd9")
    args.destination.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    patches = sorted((here / "followups").glob("*.patch"))
    if args.layout:
        patches.extend((here / "layout.patch", here / "index-pitch.patch"))
    for patch in patches:
        destination = args.destination / patch.stem
        destination.mkdir()
        for source in (header, args.cpu / "dp.cpp", args.cpu / "dp_internal.hpp"):
            shutil.copyfile(source, destination / source.name)
        subprocess.run(["patch", "-s", "-p0", "-i", str(patch)], check=True, cwd=destination)


if __name__ == "__main__":
    main()
