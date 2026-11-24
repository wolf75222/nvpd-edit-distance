"""Recreate the five historical inputs and deterministic boundary controls."""
import csv
import hashlib
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from alternatives import historical_cases

def prepare(dest):
    dest.mkdir(parents=True, exist_ok=True)
    historical = {r["case"]: r for r in csv.DictReader(Path("results/alternatives/astar.csv").open())}
    rows = []
    for name, original, (m, n) in historical_cases(dest):
        target = dest / (name + ".fa")
        if target != original:
            shutil.copyfile(original, target)
        sha = hashlib.sha256(target.read_bytes()).hexdigest()
        if sha != historical[name]["sha256"]:
            raise ValueError(f"historical input hash changed for {name}")
        rows.append((name, m, n, historical[name]["distance"], sha))
    left, right = "A"*65 + "C"*64, "C"*64 + "A"*65
    for name, a, b, distance in [
            ("periodic512", (left*4)[:512], (right*4)[:512], 128),
            ("periodic256-tail", (left*2)[:256]+"G"*256, (right*2)[:256]+"G"*256, 128),
            ("identical512", (left*4)[:512], (left*4)[:512], 0)]:
        path = dest / (name + ".fa")
        path.write_text(f">a\n{a}\n>b\n{b}\n")
        rows.append((name, len(a), len(b), distance, hashlib.sha256(path.read_bytes()).hexdigest()))
    with (dest / "cases.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow("case m n exact sha256".split())
        w.writerows(rows)

if __name__ == "__main__":
    prepare(Path(sys.argv[1]))
