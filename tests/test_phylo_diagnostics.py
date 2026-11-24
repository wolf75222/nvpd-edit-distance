import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


class PhylogenyDiagnostics(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.out = self.tmp / "output"

    def run_script(self, name, *args):
        result = subprocess.run([sys.executable, SCRIPTS / name, *args], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def rows(self, name):
        with (self.out / name).open() as stream:
            return list(csv.DictReader(stream))

    def matrix(self, labels, d):
        path = self.tmp / "matrix.tsv"
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["id", *labels])
            writer.writerows([name, *row] for name, row in zip(labels, d))
        return path

    def test_additive_unequal_branches_and_duplicate(self):
        d = [[0, 0, 2, 3, 9], [0, 0, 2, 3, 9], [2, 2, 0, 3, 9], [3, 3, 3, 0, 8], [9, 9, 9, 8, 0]]
        matrix = self.matrix(["A", "A2", "B", "C", "D"], d)
        self.run_script("phylo_fit.py", matrix, "--out", self.out)
        counts = {(r["sample"], r["criterion"]): r for r in self.rows("duplicate_sensitivity.csv")}
        self.assertEqual(counts["all", "four_point"]["violations"], "0")
        self.assertEqual(counts["distinct", "four_point"]["combinations"], "1")
        self.assertEqual(counts["distinct", "ultrametric"]["violations"], "2")
        fits = {r["tree"]: float(r["rmse"]) for r in self.rows("fit.csv")}
        self.assertAlmostEqual(fits["nj"], 0)
        self.assertGreater(fits["upgma"], 0)

    def test_rooted_clade_vs_complement_of_unrooted_split(self):
        matrix = self.matrix("ABCD", [[0, 2, 3, 9], [2, 0, 3, 9], [3, 3, 0, 8], [9, 9, 8, 0]])
        manifest = self.tmp / "groups.tsv"
        manifest.write_text("id\tgroup\nA\tx\nB\tx\nC\ty\nD\ty\n")
        self.run_script("phylo_groups.py", matrix, manifest, "--out", self.out)
        reports = {(r["tree"], r["group"]): r for r in self.rows("group_clades.csv")}
        self.assertEqual(reports["upgma", "x"]["positive_edge_match"], "1")
        self.assertEqual(reports["upgma", "y"]["positive_edge_match"], "0")
        self.assertEqual(reports["upgma", "y"]["cover_size"], "4")
        self.assertEqual(reports["nj", "y"]["positive_edge_match"], "1")

    def test_pruning_adds_lengths_and_matches_whole_fields(self):
        source = self.tmp / "tree.nwk.gz"
        with gzip.open(source, "wt") as stream:
            stream.write("(((tag|A|sample:1,tag|AA|sample:9):3,B-original:2):5,C:4):0;")
        manifest = self.tmp / "groups.tsv"
        manifest.write_text("id\tgroup\nA\tg\nB\tg\nC\th\n")
        output = self.tmp / "pruned.nwk"
        self.run_script("restrict_newick.py", source, manifest, "--alias", "B-original=B", "--out", output)
        self.assertEqual(output.read_text(), "((A:4,B:2):5,C:4);\n")
