"""Compare full-row and strided row-group OpenMP scheduling on a fixed matrix."""
import argparse
import csv
import hashlib
import subprocess
import tempfile
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=Path)
    parser.add_argument("groups", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    expected = args.matrix.read_bytes()
    matrix_sha = hashlib.sha256(expected).hexdigest()
    if matrix_sha != "fbd30e08c38b4ac27fd2d3222ab374b230fb6a3f49f6bd1bbba0a3626c252076":
        raise ValueError("expected the published 64-accession matrix")
    algorithms = ["edlib", "wfa2", "biwfa"]
    layouts = [("rows", args.rows), ("groups", args.groups)]
    with tempfile.TemporaryDirectory() as tmp, args.output.open("x", newline="") as output:
        rss = Path(tmp) / "peak-rss"
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["block", "order", "algorithm", "layout", "threads", "elapsed_s", "peak_rss_kib", "matrix_sha256"])
        for block in range(6):
            first = block % len(algorithms)
            for alg in algorithms[first:] + algorithms[:first]:
                ordered = layouts if (block + algorithms.index(alg)) % 2 == 0 else reversed(layouts)
                for order, (layout, binary) in enumerate(ordered, 1):
                    started = time.perf_counter()
                    observed = subprocess.check_output(["/usr/bin/time", "-f", "%M", "-o", rss,
                                                        binary, "matrix", "--algorithm", alg,
                                                        "--threads", "32", "--input", args.fasta])
                    elapsed = time.perf_counter() - started
                    if observed != expected:
                        raise ValueError(f"{layout}/{alg} differs from the published matrix")
                    if block:
                        writer.writerow([block, order, alg, layout, 32, f"{elapsed:.9f}", int(rss.read_text()), matrix_sha])
                        output.flush()


if __name__ == "__main__":
    main()
