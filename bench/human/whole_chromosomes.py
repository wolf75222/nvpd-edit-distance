"""Measure whole-chromosome alignments in independent, resource-limited Linux processes."""
import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path
import resource
import subprocess
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from fasta_records import read_pair


def memory_limit():
    size = 4 * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (size, size))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("manifest", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--data", default=Path("data/human"), type=Path)
    args = p.parse_args()
    fields = "job case method trial m n status distance seconds wall_seconds peak_rss_kib input_sha256".split()
    with args.manifest.open() as manifest:
        entries = [r for r in csv.DictReader(manifest, delimiter="\t") if r["case"].endswith("-assemblies")]
    pairs = {r["case"]: r for r in entries}
    if not pairs or len(pairs) != len(entries):
        raise ValueError("expected unique assembly pairs in the manifest")
    with args.out.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for case, entry in pairs.items():
            fasta = args.data / (case + ".fa")
            digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
            if digest != entry["sha256"]:
                raise ValueError(case)
            expected_lengths = int(entry["m"]), int(entry["n"])
            if tuple(map(len, read_pair(fasta))) != expected_lengths:
                raise ValueError(f"FASTA dimensions disagree with manifest for {case}")
            scores = set()
            for trial in range(1, 4 if case == "chr20-assemblies" else 2):
                for method in (("edlib", "wfa2") if trial % 2 else ("wfa2", "edlib")):
                    row = dict(job=os.environ["SLURM_JOB_ID"], case=case, method=method, trial=trial,
                               m=entry["m"], n=entry["n"], input_sha256=digest)
                    start = time.monotonic()
                    try:
                        proc = subprocess.run([args.binary,"pair","--input",fasta,"--algorithm",method],
                                              capture_output=True, text=True, timeout=120, preexec_fn=memory_limit)
                        if proc.returncode:
                            row["status"] = f"exit_{proc.returncode}"
                            print(case, method, proc.stderr[-400:], flush=True)
                        else:
                            value = json.loads(proc.stdout)
                            if (value["m"], value["n"]) != expected_lengths:
                                raise ValueError(f"FASTA dimensions disagree with manifest for {case}")
                            row.update(status="ok", **{k:value[k] for k in ("distance","seconds","peak_rss_kib")})
                            scores.add(value["distance"])
                    except subprocess.TimeoutExpired:
                        row["status"] = "wall_limit_120s"
                    row["wall_seconds"] = time.monotonic() - start
                    writer.writerow(row)
                    stream.flush()
            if len(scores) > 1:
                raise ValueError(f"discordant distances on {case}")


if __name__ == "__main__":
    main()
