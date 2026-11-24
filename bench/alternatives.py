"""Run exact single-core comparisons from fresh processes (project root)."""
import argparse
import csv
import hashlib
import json
import random
import subprocess
import tempfile
from pathlib import Path


def fasta_pair(path):
    records = path.read_text().split(">")[1:]
    return ["".join(record.splitlines()[1:]) for record in records]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--machine", required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    variants = ["banded", "banded-trim", "quicked", "wfa2", "edlib"]
    fields = "case variant trial m n distance seconds peak_rss_kib machine sha256".split()
    order = random.Random(2026)
    with tempfile.TemporaryDirectory() as tmp, args.out.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        cases = [(f"20k-{rate}", next(Path("data/synthetic/pairs").glob(f"n20000_r{rate}_p1_*.fa")))
                 for rate in ["0.02", "0.05", "0.1"]]
        a, b = fasta_pair(Path("data/human/chr21-2m.fa"))
        def pair(name, x, y):
            path = Path(tmp, name + ".fa")
            path.write_text(f">a\n{x}\n>b\n{y}\n")
            cases.append((name, path))
        pair("chr21-50k", a[:50000], b[:50000])
        # Human reference segment with a dispersed-within-block substitution control.
        a = a[:100000]
        core = a[49936:50064].translate(str.maketrans("ACGT", "CGTA"))
        pair("chr21-center128", a, a[:49936] + core + a[50064:])
        for case, path in cases:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            scores = set()
            for trial in range(1, 6):
                order.shuffle(variants)
                for variant in variants:
                    cmd = ([args.build / "editdist-quicked", path] if variant == "quicked" else
                           [args.build / "editdist", "pair", "--input", path, "--algorithm", variant])
                    result = json.loads(subprocess.check_output(cmd, timeout=180))
                    scores.add(result["distance"])
                    if len(scores) != 1:
                        raise ValueError(f"distance disagreement: {case} {variant} {scores}")
                    writer.writerow({k: result[k] for k in fields if k in result} |
                                    dict(case=case, variant=variant, trial=trial,
                                         machine=args.machine, sha256=digest))
                    stream.flush()


if __name__ == "__main__":
    main()
