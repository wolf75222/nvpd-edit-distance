"""Compare full-matrix and linear-space traceback on prefixes of one synthetic pair."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess
import tempfile


def replay(a, b, cigar, score):
    i = j = cost = end = 0
    for run in re.finditer(r"([1-9][0-9]*)([=XID])", cigar):
        if run.start() != end:
            raise ValueError("invalid CIGAR")
        end = run.end()
        count, op = int(run[1]), run[2]
        if op in "=X":
            first, second = a[i:i+count], b[j:j+count]
            if len(first) != count or len(second) != count:
                raise ValueError("CIGAR exceeds sequence length")
            if any((x == y) != (op == "=") for x, y in zip(first, second)):
                raise ValueError("CIGAR symbol mismatch")
        i += count if op != "I" else 0
        j += count if op != "D" else 0
        cost += count if op != "=" else 0
    if (end, i, j, cost) != (len(cigar), len(a), len(b), score):
        raise ValueError("CIGAR length or cost mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("pair", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    records = [b"".join(record.splitlines()[1:]) for record in args.pair.read_bytes().split(b">")[1:]]
    if len(records) != 2 or min(map(len, records)) < 25000:
        raise ValueError("expected a pair of at least 25000 bases")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(2026)
    with tempfile.TemporaryDirectory() as tmp, args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, "case algorithm trial m n distance seconds peak_rss_kib sha256".split(), lineterminator="\n")
        writer.writeheader()
        path = Path(tmp) / "prefix.fa"
        for length in (2000, 8000, 25000):
            a, b = (record[:length] for record in records)
            fasta = b">A\n" + a + b"\n>B\n" + b + b"\n"
            path.write_bytes(fasta)
            argv = [args.binary.resolve(), "pair", "--input", path, "--threads", "1", "--algorithm"]
            oracle = json.loads(subprocess.check_output([*argv, "edlib"]))["distance"]
            for trial in range(1, 6):
                methods = ["hirschberg", "seq"] if length <= 20000 else ["hirschberg"]
                rng.shuffle(methods)
                for method in methods:
                    row = json.loads(subprocess.check_output([*argv, method, "--trace"]))
                    replay(a, b, row["cigar"], oracle)
                    if row["distance"] != oracle or (row["m"], row["n"]) != (length, length):
                        raise ValueError("score or input dimensions differ from reference")
                    writer.writerow(dict(case=f"prefix-{length}", algorithm=method, trial=trial,
                                         sha256=hashlib.sha256(fasta).hexdigest(),
                                         **{k: row[k] for k in ("m","n","distance","seconds","peak_rss_kib")}))
                    stream.flush()


if __name__ == "__main__":
    main()
