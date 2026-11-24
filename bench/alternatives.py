"""Run exact single-core comparisons from fresh processes (project root)."""
import argparse
import csv
import hashlib
import json
import random
import subprocess
import tempfile
from pathlib import Path

from fasta_records import read_pair


def fasta_pair(path):
    return [s.decode("ascii") for s in read_pair(path)]


def historical_cases(tmp):
    # The generator can leave older seeds in pairs/. Only the manifest
    # identifies the active inputs for this campaign.
    manifest = Path("data/synthetic/manifest.tsv")
    rates = ("0.02", "0.05", "0.1")
    with manifest.open(newline="") as f:
        selected = [r for r in csv.DictReader(f, delimiter="\t")
                    if r["length"] == "20000" and r["rep"] == "1" and r["rate"] in rates]
    inputs = {r["rate"]: r for r in selected}
    if len(selected) != len(rates) or set(inputs) != set(rates):
        raise ValueError("expected one 20k replicate per error rate in the synthetic manifest")
    cases = [(f"20k-{rate}", manifest.parent / inputs[rate]["path"],
              (20000, int(inputs[rate]["mutant_length"]))) for rate in rates]
    a, b = fasta_pair(Path("data/human/chr21-2m.fa"))
    if len(a) < 100000 or len(b) < 50000:
        raise ValueError("human FASTA is too short for the specified windows")

    def pair(name, x, y):
        path = Path(tmp, name + ".fa")
        path.write_text(f">a\n{x}\n>b\n{y}\n")
        cases.append((name, path, (len(x), len(y))))
    pair("chr21-50k", a[:50000], b[:50000])
    # Human reference segment with a dispersed-within-block substitution control.
    a = a[:100000]
    core = a[49936:50064].translate(str.maketrans("ACGT", "CGTA"))
    pair("chr21-center128", a, a[:49936] + core + a[50064:])
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--machine", required=True)
    parser.add_argument("--astar-binary", type=Path, help="benchmark both modes of the upstream A*PA2 library")
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    variants = ["banded", "banded-trim", "quicked", "wfa2", "edlib"]
    if args.astar_binary:
        variants.extend(["astarpa2-simple", "astarpa2-full"])
    fields = "case variant trial m n distance seconds peak_rss_kib machine sha256".split()
    order = random.Random(2026)
    with tempfile.TemporaryDirectory(dir=args.out.parent) as tmp:
        pending = Path(tmp, "trials.csv")
        with pending.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            cases = historical_cases(tmp)
            for case, path, lengths in cases:
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                if tuple(map(len, read_pair(path))) != lengths:
                    raise ValueError(f"unexpected sequence lengths in {case}")
                scores = set()
                for trial in range(1, 6):
                    order.shuffle(variants)
                    for variant in variants:
                        cmd = ([args.build / "editdist-quicked", path] if variant == "quicked" else
                               [args.build / "editdist", "pair", "--input", path, "--algorithm", variant])
                        if variant.startswith("astarpa2-"):
                            cmd = [args.astar_binary, variant.split("-")[-1], path]
                        result = json.loads(subprocess.check_output(cmd, timeout=180))
                        if (result["m"], result["n"]) != lengths:
                            raise ValueError(f"FASTA dimensions disagree for {case}, {variant}")
                        scores.add(result["distance"])
                        if len(scores) != 1:
                            raise ValueError(f"distance disagreement: {case} {variant} {scores}")
                        writer.writerow({k: result[k] for k in fields if k in result} |
                                        dict(case=case, variant=variant, trial=trial,
                                             machine=args.machine, sha256=digest))
                        stream.flush()
        pending.replace(args.out)


if __name__ == "__main__":
    main()
