"""Data rows for the CPU, accelerator and distributed comparison tables."""

import csv
from collections import defaultdict
from statistics import median


def grouped(root, path, keys, keep=lambda row: True):
    groups = defaultdict(list)
    with (root / path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            if keep(row):
                groups[tuple(row[key] for key in keys)].append(row)
    return groups


def time(rows):
    return median(float(row["seconds"]) for row in rows)


def profile_tables(root, number, table_row, macro):
    result = []

    def table(name, rows):
        result.append(macro(name, "\n".join(table_row(*row) for row in rows)))

    ab = grouped(
        root, "profile/cpu/ab.csv", ("case", "threads", "variant"), lambda r: r["comparison"] == "baseline-v1"
    )
    rows = []
    for shape in ("square", "wide", "over_l3"):
        values = []
        for p in ("1", "8", "16"):
            orig = {r["trial"]: r for r in ab[shape, p, "baseline"]}
            simd = ab[shape, p, "v1"]
            assert len(orig) == len(simd) == 3
            assert all(orig[r["trial"]]["input_fnv64"] == r["input_fnv64"] for r in simd)
            values.append(
                number(median(float(orig[r["trial"]]["seconds"]) / float(r["seconds"]) for r in simd))
            )
        case = simd[0]
        rows.append(
            [r"\(" + number(int(case["m"]), 0) + r"\times" + number(int(case["n"]), 0) + r"\)", *values]
        )
    table("CpuPairedRows", rows)

    scaling = grouped(root, "profile/cpu/scaling.csv", ("study", "threads"), lambda r: r["variant"] == "v1")
    rows = []
    for p in ("1", "2", "4", "8", "16"):
        strong, weak = (time(scaling[study, p]) for study in ("strong", "weak"))
        speed = time(scaling["strong", "1"]) / strong
        rows.append(
            [
                p,
                number(1000 * strong),
                number(speed),
                number(100 * speed / int(p), 1),
                number(1000 * weak),
                number(100 * time(scaling["weak", "1"]) / weak, 1),
            ]
        )
    table("CpuScalingRows", rows)

    rows = []
    for rate, job in ((2, "838170"), (10, "838173")):
        gpu = grouped(
            root, f"profile/gpu/{job}/trials.csv", ("m", "variant", "run"), lambda r: r["mode"] == "warm"
        )
        row = [str(rate) + r"\,\%"]
        for n in (5000, 20000, 100000, 250000, 500000, 1000000):
            if (str(n), "baseline", "1") not in gpu:
                row.append("--")
                continue
            baseline, candidate = (
                median(time(gpu[str(n), v, str(i)]) for i in (1, 2, 3)) for v in ("baseline", "candidate")
            )
            row.append(number(baseline / candidate))
        rows.append(row)
    table("GpuAdaptiveRows", rows)

    alternatives = grouped(root, "alternatives/trials.csv", ("case", "variant"))
    rows = []
    for case, label in (
        ("20k-0.02", r"20\,000, 2\,\%"),
        ("20k-0.05", r"20\,000, 5\,\%"),
        ("20k-0.1", r"20\,000, 10\,\%"),
        ("chr21-50k", r"chr21, 50\,000"),
        ("chr21-center128", "chr21, contrôle"),
    ):
        score = int(alternatives[case, "edlib"][0]["distance"])
        row = [label, number(score, 0)]
        for method in ("banded", "banded-trim", "edlib", "wfa2", "quicked"):
            sample = alternatives[case, method]
            assert len(sample) == 5 and all(int(r["distance"]) == score for r in sample)
            row.append(number(1000 * time(sample), 3))
        rows.append(row)
    table("AlternateRows", rows)

    mpi = grouped(root, "mpi/838226/trials.csv", ("case", "nodes"))
    rows = []
    for case, label in (
        ("square", "Paire quasi carrée"),
        ("wide", "Paire large"),
        ("synthetic64", r"Matrice, 64 chaînes de 12\,000 bases"),
        ("phix64", "Matrice, 64 génomes de phage"),
    ):
        single, double = (time(mpi[case, n]) for n in ("1", "2"))
        rows.append(
            [
                label,
                number(single, 4),
                number(double, 4),
                number(single / double, 3),
                number(50 * single / double, 1),
            ]
        )
    table("MpiRows", rows)

    human = grouped(root, "human/838220/trials.csv", ("case", "algorithm"), lambda r: r["status"] == "ok")
    rows = []
    for case, label in (
        ("chr21-2m", "Fenêtre 2 M"),
        ("chr21-4m", "Fenêtre 4 M"),
        ("chr21-full-insert64", "Contrôle entier"),
    ):
        ref = human[case, "edlib"][0]
        dims = r"\(" + number(int(ref["m"]), 0) + r"\times" + number(int(ref["n"]), 0) + r"\)"
        row = [label, dims, number(int(ref["distance"]), 0)]
        for method in ("edlib", "wfa2", "nvpd-cuda"):
            sample = human.get((case, method))
            row.append(number(time(sample), 3) if sample else "--")
        rows.append(row)
    table("HumanRows", rows)
    return result
