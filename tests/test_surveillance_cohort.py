"""Check agreement of the locked cohort, published distances and report appendix."""
import csv
import hashlib
import re
import unittest
from collections import Counter
from itertools import combinations
from pathlib import Path

from scripts.datasets.ncbi_metadata import LINEAGES, REFERENCE
from scripts.phylogeny.tree import read_matrix

ROOT = Path(__file__).resolve().parents[1]


def rows(filename):
    with (ROOT / "results/ncbi" / filename).open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


class SurveillanceCohort(unittest.TestCase):
    def test_accessions_and_count_bound(self):
        members = rows("sars_metadata.tsv")
        content = rows("sars_content.tsv")
        ids, matrix = read_matrix(ROOT / "results/phylogeny/sars-edit/matrix.tsv")
        self.assertEqual(ids, [m["id"] for m in members])
        self.assertEqual(ids, [m["id"] for m in content])
        self.assertEqual(len(set(ids)), 64)
        self.assertEqual([m["id"] for m in members if m["group"] == "reference"], [REFERENCE])
        self.assertEqual(Counter(m["group"] for m in members), {"reference": 1, **dict.fromkeys(LINEAGES, 9)})
        for m in members:
            self.assertRegex(m["id"], r"^[A-Z_]+\d+\.\d+$")
            if m["group"] != "reference":
                self.assertEqual(m["group"], m["pango"])

        appendix = (ROOT / "report/sections/accessions.tex").read_text().replace(r"\_", "_")
        printed = re.findall(r"\\texttt\{([A-Z_]+\d+\.\d+)\} & ([\d\\,]+)", appendix)
        self.assertEqual([(key, count.replace(r"\,", "")) for key, count in printed],
                         [(m["id"], m["n"]) for m in content])
        n = [int(m["n"]) for m in content]
        self.assertEqual(sum(n), 65565)
        self.assertEqual(sum(int(m["ambiguous"]) for m in content), 65588)
        bounds = [(abs(n[i] - n[j]), matrix[i][j]) for i, j in combinations(range(len(ids)), 2)]
        self.assertTrue(all(bound <= distance for bound, distance in bounds))
        self.assertEqual(sum(2 * bound >= distance for bound, distance in bounds), 915)

    def test_benchmark_matrix_digests(self):
        base = ROOT / "results/phylogeny/sars-edit"
        fingerprint = hashlib.sha256((base / "matrix.tsv").read_bytes()).hexdigest()
        with (base / "verification.csv").open() as source:
            for run in csv.DictReader(source):
                self.assertEqual((run["records"], run["pairs"], run["matrix_sha256"]),
                                 ("64", "2016", fingerprint))
