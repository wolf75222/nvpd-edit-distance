"""Cross-check the published score list with WFA2-lib on the same window files."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def verify(inputs, scores, binary):
    with scores.open() as source:
        expected = list(csv.DictReader(source,delimiter="\t"))
    checked = 0
    with tempfile.TemporaryDirectory() as directory:
        pair = Path(directory) / "window.fa"
        for name in ("matched-768","pair-64k","pair-128k"):
            lines = (inputs / (name+".fa")).read_bytes().splitlines()
            for i in range(0,len(lines),4):
                a, x, b, y = lines[i:i+4]
                row = expected[checked]
                if (a[1:].decode(),b[1:].decode()) != (row["id_a"],row["id_b"]):
                    raise ValueError("window identifiers disagree")
                pair.write_bytes(b"\n".join((a,x,b,y))+b"\n")
                result = subprocess.run([binary,"pair","--input",pair,"--algorithm","wfa2"],
                                        capture_output=True,check=True,text=True,timeout=120)
                if json.loads(result.stdout)["distance"] != int(row["edlib"]):
                    raise ValueError(f"edlib and WFA2 disagree at {row['id_a']}")
                checked += 1
    if checked != len(expected):
        raise ValueError("unverified scores remain")
    print(f"{checked} WFA2 scores match the edlib reference")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("inputs","scores","binary"): parser.add_argument(field,type=Path)
    opts = parser.parse_args()
    verify(opts.inputs,opts.scores,opts.binary.resolve())
