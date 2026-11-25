"""Exercise the CPU routing policies through their public FASTA entry point."""

from collections import Counter
import csv
import io
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def hash_words(sequence):
    """Hash every distinct 15-mer for comparison with the streaming sketch."""
    hashes = set()
    mask = (1 << 64) - 1
    digits = sequence.translate(str.maketrans("ACGT", "0123"))
    for start in range(len(digits) - 14):
        x = int(digits[start:start + 15], 4) + 0x9e3779b97f4a7c15
        x = ((x ^ (x >> 30)) * 0xbf58476d1ce4e5b9) & mask
        x = ((x ^ (x >> 27)) * 0x94d049bb133111eb) & mask
        hashes.add(x ^ (x >> 31))
    return hashes


def run(binary, work, records, policy, threads):
    fasta, matrix, decisions = (work / x for x in ("input.fa", "matrix.tsv", "decisions.tsv"))
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    output = subprocess.check_output(
        [binary, "--policy", policy, "--input", fasta, "--threads", str(threads),
         "--warmups", "1", "--repeats", "2", "--matrix", matrix, "--decisions", decisions],
        text=True,
    )
    timings = list(csv.DictReader(io.StringIO(output)))
    labels = [name for name, _ in records]
    n = len(labels)
    pairs = n * (n - 1) // 2
    assert [row["phase"] for row in timings] == ["warmup", "measure", "measure"]
    for row in timings:
        assert row["policy"] == policy and int(row["threads"]) == threads
        assert int(row["pairs"]) == int(row["cpu_pairs"]) == pairs
        assert int(row["gpu_pairs"]) == 0
        elapsed, estimate = float(row["seconds"]), float(row["sketch_seconds"])
        assert math.isfinite(elapsed) and 0 <= estimate <= elapsed

    with matrix.open() as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    assert rows[0][1:] == labels and [r[0] for r in rows[1:]] == labels
    scores = {}
    for i, row in enumerate(rows[1:]):
        assert len(row) == n + 1 and int(row[i + 1]) == 0
        for j in range(i):
            score = int(row[j + 1])
            assert score >= 0 and score == int(rows[j + 1][i + 1])
            scores[tuple(sorted((labels[i], labels[j])))] = score
    policy_decisions = {}
    if policy == "cpu-sketch":
        with decisions.open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                key = tuple(sorted((row["id_i"], row["id_j"])))
                assert key not in policy_decisions and int(row["distance"]) == scores[key]
                policy_decisions[key] = (row["class"], int(row["shared"]), int(row["sampled"]))
        assert policy_decisions.keys() == scores.keys()
    return scores, policy_decisions


def check_composition_routing(binary, work):
    # These pairs share all 15-mers; N, R and Y differ only in the ambiguous suffix.
    stem = "A" * 80
    inputs = dict(ref=stem + "N" * 20, below=stem + "N" * 16 + "R" * 4,
                  at=stem + "N" * 15 + "R" * 5, above=stem + "N" * 14 + "Y" * 6,
                  short="NR")
    entries = list(inputs.items())
    oracle, _ = run(binary, work, entries, "cpu-edlib", 1)
    scores, _ = run(binary, work, entries, "hybrid-counts", 4)
    assert scores == oracle
    with (work / "decisions.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            a, b = inputs[row["id_i"]], inputs[row["id_j"]]
            x, y = Counter(a), Counter(b)
            bound = max(sum((x - y).values()), sum((y - x).values()))
            assert int(row["lower_bound"]) == bound <= int(row["distance"])
            assert row["class"] == ("fallback" if "NR" in (a, b) else "close")
            expected = "edlib" if row["class"] == "fallback" or 20 * bound >= max(len(a), len(b)) else "wfa2"
            assert row["cpu_library_choice"] == expected


def main(binary):
    rng = random.Random(8341)
    ref = "".join(rng.choices("ACGT", k=1024))
    unrelated = "".join(rng.choices("ACGT", k=1024))
    mosaic = ref[:700] + unrelated[700:]
    full_ref, full_mosaic = hash_words(ref), hash_words(mosaic)
    sample = set(sorted(full_ref | full_mosaic)[:256])
    shared = len(sample & full_ref & full_mosaic)
    assert 64 < shared < 256
    four_words = "N".join(base * 15 for base in "ACGT")
    records = [("ref", ref), ("copy", ref), ("mosaic", mosaic), ("other", unrelated), ("A20", "A" * 20),
               ("A21", "A" * 21), ("C20", "C" * 20), ("short", "AC"),
               ("quart", four_words), ("fifth", four_words + "N" + "AC" * 7 + "G"),
               ("break", "A" * 8 + "N" + "A" * 8), ("R", "R"), ("N", "N")]
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        check_composition_routing(binary, work)
        oracle, _ = run(binary, work, records, "cpu-edlib", 1)
        assert oracle[("copy", "ref")] == 0 and oracle[("N", "R")] == 1
        reference_routes = None
        for order in (records, records[::-1]):
            for threads in (1, 4):
                for policy in ("cpu-edlib", "cpu-wfa2", "cpu-sketch"):
                    scores, routes = run(binary, work, order, policy, threads)
                    assert scores == oracle, (policy, threads)
                    if not routes:
                        continue
                    assert routes[("copy", "ref")] == ("close", 256, 256)
                    assert routes[("mosaic", "ref")] == ("close", shared, 256)
                    assert routes[("other", "ref")][0] == "far"
                    assert routes[("A20", "A21")] == ("close", 1, 1)
                    assert routes[("A20", "C20")] == ("far", 0, 2)
                    assert routes[("A20", "quart")] == ("close", 1, 4)
                    assert routes[("A20", "fifth")] == ("far", 1, 5)
                    for key, (decision, _, _) in routes.items():
                        if {"short", "break", "N", "R"}.intersection(key):
                            assert decision == "fallback", key
                    if reference_routes is None:
                        reference_routes = routes
                    assert routes == reference_routes


if __name__ == "__main__":
    main(sys.argv[1])
