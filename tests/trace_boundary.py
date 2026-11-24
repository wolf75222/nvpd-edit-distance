"""Exercise the dense traceback limit on square inputs through the public CLI."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bench.cpu.traceback.traceback import replay
from scripts.datasets.synthetic import make_pair, pair_seed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("backends", nargs="+", choices=("seq", "nvpd-omp", "nvpd-cuda",
                                                         "nvpd-cuda-local", "nvpd-kokkos"))
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    seed = pair_seed(2026, 100_000, "0.02", 1)
    reference, mutant, _ = make_pair(100_000, 2_000, seed)
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(("m", "n", "input_sha256", "algorithm", "mode", "outcome", "distance",
                     "device_bytes"))
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "pair.fa"
        cli = [args.binary.resolve(), "pair", "--input", path, "--threads", str(args.threads),
               "--algorithm"]
        for m, n in ((20_000, 20_000), (20_001, 20_000), (20_000, 20_001)):
            a, b = reference[:m], mutant[:n]
            data = b">A\n" + a + b"\n>B\n" + b + b"\n"
            path.write_bytes(data)
            digest = hashlib.sha256(data).hexdigest()
            oracle = json.loads(subprocess.check_output([*cli, "edlib"]))["distance"]
            backends = [(method, (False, True)) for method in args.backends]
            backends.extend([("diagonal", (False,)), ("hirschberg", (True,))])
            for backend, modes in backends:
                for trace in modes:
                    mode = "trace" if trace else "score"
                    cmd = [*cli, backend, *(["--trace"] if trace else [])]
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    if trace and backend in args.backends and max(m, n) > 20_000:
                        if result.returncode == 0 or "traceback supports at most 20000" not in result.stderr:
                            raise AssertionError((backend, m, n, result.stderr))
                        writer.writerow((m, n, digest, backend, mode, "length-limit", "", ""))
                        continue
                    result.check_returncode()
                    row = json.loads(result.stdout)
                    if (row["m"], row["n"], row["distance"]) != (m, n, oracle):
                        raise AssertionError((backend, m, n, row, oracle))
                    if trace:
                        replay(a, b, row["cigar"], oracle)
                    writer.writerow((m, n, digest, backend, mode, "pass", oracle, row["device_bytes"]))


if __name__ == "__main__":
    main()
