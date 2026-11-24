"""Every wall-clock sample in the published surveillance comparison covers the same matrix."""
import csv
import hashlib
from collections import defaultdict
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phylogeny"


class SurveillanceTimings(unittest.TestCase):
    def test_complete_blocks_and_reference(self):
        fingerprint = hashlib.sha256((DATA / "sars-edit/matrix.tsv").read_bytes()).hexdigest()
        expected = {
            "839288": {(name, width) for name in ("edlib", "wfa2", "nvpd-split", "nvpd-adaptive")
                       for width in (1, 32)} | {("nvpd-omp", 32)},
            "839297": {(name, width) for name in ("edlib", "wfa2", "biwfa", "nvpd-split")
                       for width in (1, 32)},
            "839289": {("edlib", 1), ("wfa2", 1), ("nvpd-cuda-local", 1)}
        }
        expected["839304"] = expected["839297"]
        for job, configs in expected.items():
            blocks = defaultdict(list)
            with (DATA / "sars-romeo" / job / "trials.csv").open() as source:
                for row in csv.DictReader(source):
                    self.assertEqual(row["job"], job)
                    self.assertEqual(row["matrix_sha256"], fingerprint)
                    self.assertGreater(float(row["wall_s"]), 0)
                    blocks[int(row["block"])].append((row["algorithm"], int(row["threads"])))
            self.assertEqual(set(blocks), set(range(6)), job)
            for cases in blocks.values():
                self.assertCountEqual(cases, configs)

    def test_scheduling_control(self):
        fingerprint = hashlib.sha256((DATA / "sars-edit/matrix.tsv").read_bytes()).hexdigest()
        blocks = defaultdict(list)
        with (DATA / "sars-schedule/839293/trials.csv").open() as source:
            for row in csv.DictReader(source):
                self.assertEqual(row["matrix_sha256"], fingerprint)
                blocks[int(row["block"])].append((row["algorithm"], row["layout"]))
        self.assertEqual(set(blocks), set(range(1, 6)))
        for cases in blocks.values():
            self.assertCountEqual(cases, [(alg, layout) for alg in ("edlib", "wfa2", "biwfa")
                                         for layout in ("rows", "groups")])

    def test_literal_scores_in_the_gpu_dispatcher(self):
        with (ROOT / "results/ncbi/sars_content.tsv").open() as source:
            unambiguous = sum(int(r["ambiguous"]) == 0 for r in csv.DictReader(source, delimiter="\t"))
        eligible = unambiguous * (unambiguous - 1) // 2
        with (ROOT / "results/wfa-gpu/839298/trials.csv").open() as source:
            blocks = defaultdict(list)
            for row in csv.DictReader(source):
                blocks[int(row["block"])].append(row["method"])
                self.assertEqual(int(row["wrong_scores"]), 0)
                self.assertEqual(int(row["pairs"]), 2016)
                self.assertEqual(int(row["cpu_scored"]) + int(row["gpu_scored"]), 2016)
                if row["method"] == "hybrid":
                    self.assertEqual(int(row["gpu_submitted"]), eligible)
                    self.assertEqual(int(row["gpu_scored"]), eligible)
        self.assertEqual(set(blocks), set(range(6)))
        for cases in blocks.values():
            self.assertCountEqual(cases, ["hybrid", "wfa2", "edlib"])
