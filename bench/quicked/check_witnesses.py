"""Record raw and certified scores for the three fixed controls on another ISA."""
import csv
from pathlib import Path
import subprocess
import sys

build, inputs, machine, output = sys.argv[1:]
with Path(output).open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow("platform case variant score exact correct".split())
    for case in ("periodic512", "periodic256-tail", "identical512"):
        path = Path(inputs) / (case + ".fa")
        def score(exe, variant):
            return int(subprocess.check_output([Path(build)/exe,path,variant,"1"], text=True).split(",")[1])
        exact = score("measure-quicked-native", "edlib")
        for exe, variant in (("measure-quicked", "quicked-certified"), ("measure-quicked-stock", "quicked-stock")):
            value = score(exe, variant)
            if variant == "quicked-certified" and value != exact:
                raise ValueError(case)
            writer.writerow((machine, case, variant, value, exact, int(value==exact)))
