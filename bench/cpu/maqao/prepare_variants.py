"""Materialize the five-way EPYC experiment from the original NvPD source."""
import argparse
import hashlib
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path, help="checkout of ba97a3d")
    parser.add_argument("destination", type=Path, help="new directory for candidate sources")
    args = parser.parse_args()
    original = (args.baseline / "src/cpu/dp.cpp").read_bytes()
    expected = "8eb5500a8977a3ed5da809189a8c77e2caf410ec9fc7840e18ec2952519947c7"
    if hashlib.sha256(original).hexdigest() != expected:
        raise ValueError("unexpected baseline source; use the checkout of ba97a3d")
    args.destination.mkdir(parents=True)
    for name in ("centered.cpp", "blocked.cpp"):
        (args.destination / name).write_bytes(original)
    patch = Path(__file__).with_name("variants.patch").resolve()
    subprocess.run(["patch", "-s", "-p0", "-i", str(patch)], cwd=args.destination, check=True)


if __name__ == "__main__":
    main()
