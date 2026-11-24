#!/usr/bin/env python3
"""Refresh measurements.tex from the downloaded ROMEO CSVs."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median
from phylograms import phylograms
from profile_tables import profile_tables

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "results" / "romeo"
RATES = (0.02, 0.05, 0.10)
METHODS = (("DP séquentielle", "seq", 1), ("Anti-diagonales, 8 fils", "diagonal", 8),
           ("NvPD, 8 fils", "nvpd-omp", 8), ("NvPD, GPU", "nvpd-cuda", 1),
           ("edlib, 1 fil", "edlib", 1), ("WFA2-lib, 1 fil", "wfa2", 1))


def number(x, decimals=2):
    if x is None:
        return "--"
    return f"{x:,.{decimals}f}".replace(",", r"\,").replace(".", "{,}")


def macro(name, value):
    return r"\newcommand{\%s}{%s}" % (name, value)


def table_row(*cells):
    return " & ".join(map(str, cells)) + r" \\"


def main():
    rows = []
    for path in sorted(DATA.glob("*-*.csv")):
        with path.open() as f:
            rows.extend(r for r in csv.DictReader(f) if r["status"] == "ok")
    if not rows:
        raise ValueError(f"no successful trials in {DATA}")
    scores = defaultdict(set)
    for r in rows:
        scores[r["input_sha256"]].add(int(r["distance"]))
    if any(len(s) > 1 for s in scores.values()):
        raise ValueError("algorithms disagree for at least one synthetic pair")

    def sample(length, rate, algo, threads, job=None):
        return [r for r in rows if int(r["length"]) == length and float(r["rate"]) == rate
                and r["algorithm"] == algo and int(r["threads"]) == threads
                and (job is None or r["job"] == job)]

    # Fixed pilot jobs have three independent process launches on replicate 1.
    def pilot(algo, threads=1):
        job = "837262" if algo == "nvpd-cuda" else "837261"
        rs = sample(20000, 0.05, algo, threads, job)
        if len(rs) != 3 or {r["rep"] for r in rs} != {"1"}:
            raise ValueError(f"incomplete pilot for {algo}/{threads}")
        return median(float(r["seconds"]) for r in rs)

    lines, scaling = [], []
    for p in (1, 2, 4, 8):
        anti, nv = pilot("diagonal", p), pilot("nvpd-omp", p)
        scale = pilot("nvpd-omp") / nv
        scaling.append(table_row(p, number(1e3 * anti, 1), number(pilot("diagonal") / anti),
                                 number(1e3 * nv, 1), number(scale), number(100 * scale / p, 0)))
    lines += [macro("ScalingRows", "\n".join(scaling)), macro("RefMs", number(1e3 * pilot("seq"), 1)),
              macro("PilotGPU", number(1e3 * pilot("nvpd-cuda"), 1)),
              macro("PilotSpeed", number(pilot("seq") / pilot("nvpd-omp", 8)))]

    long_rows = []
    for name, algo, threads in METHODS:
        cells = [name]
        selected = []
        for rate in RATES:
            rs = sample(1000000, rate, algo, threads)
            selected.extend(rs)
            pairs = defaultdict(list)
            for r in rs:
                pairs[r["input_sha256"]].append(float(r["seconds"]))
            if pairs:
                t = median(median(times) for times in pairs.values())
                cells.append(number(t) + " (" + str(len(pairs)) + ")")
            else:
                cells.append("--")
        mid = [r for r in selected if float(r["rate"]) == 0.05]
        pair_rates = defaultdict(list)
        for r in mid:
            pair_rates[r["input_sha256"]].append(int(r["m"]) * int(r["n"]) / (1e9 * float(r["seconds"])))
        gcups = median(median(v) for v in pair_rates.values()) if pair_rates and algo not in ("edlib", "wfa2") else None
        rss = max((int(r["peak_rss_kib"]) / 1024 for r in selected), default=None)
        device = max((int(r["device_bytes"]) / 2**20 for r in selected if int(r["device_bytes"])), default=None)
        long_rows.append(table_row(*cells, number(gcups, 1), number(rss, 1), number(device, 1)))
    lines.append(macro("LongRows", "\n".join(long_rows)))
    seq = [r for r in sample(1000000, 0.05, "seq", 1) if r["rep"] == "1"]
    speed_note = ""
    if seq:
        digest = seq[0]["input_sha256"]
        gpu = [r for r in sample(1000000, 0.05, "nvpd-cuda", 1) if r["input_sha256"] == digest]
        if gpu:
            ratio = median(float(r["seconds"]) for r in seq) / median(float(r["seconds"]) for r in gpu)
            speed_note = (r" À 5\,\%, sur la première paire, l'accélération du GPU par rapport à la DP "
                          "séquentielle vaut " + number(ratio, 1) + ".")
    lines.append(macro("LongSpeed", speed_note))
    cohort = []
    for job in ("837261", "837271"):
        with (DATA / f"cohort-{job}" / "timings.csv").open() as f:
            cohort.extend(csv.DictReader(f))
    if len({(r["fasta_sha256"], r["matrix_sha256"]) for r in cohort}) != 1:
        raise ValueError("cohort matrix digests differ")
    for r in cohort:
        if r["job"] == "837261":
            key = {"edlib": "PhageEdlib", "wfa2": "PhageWfa"}[r["algorithm"]]
            lines.append(macro(key, number(float(r["wall_seconds"]) * 1000, 1)))
        elif r["algorithm"] == "nvpd-cuda":
            lines.append(macro("PhageGpu", number(float(r["wall_seconds"]), 1)))
    lines.extend(profile_tables(DATA.parent, number, table_row, macro))
    (HERE / "measurements.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    figure = phylograms(DATA / "cohort-837261")
    (HERE / "figures" / "arbres.tex").write_text(figure, encoding="utf-8")


if __name__ == "__main__":
    main()
