import csv
import hashlib
import io
import tempfile
import unittest
from decimal import Decimal
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

from scripts.datasets import human_reference
from scripts.datasets.synthetic import generate
from scripts.phylogeny.tree import build_tree, read_matrix


class DatasetTests(unittest.TestCase):
    def test_chromosome_cache_recovers_only_with_verified_content(self):
        reference = b">example.1\nac\ngt\n"
        record = ("example.1", hashlib.sha256(b"ACGT").hexdigest())
        with tempfile.TemporaryDirectory() as tmp, patch.dict(human_reference.REFS, {"sample": record}):
            cache = Path(tmp)
            entry = cache / "sample.fa"
            with patch("urllib.request.urlopen", return_value=io.BytesIO(b"")):
                with self.assertRaisesRegex(ValueError, "FASTA accession"):
                    human_reference.chromosome(cache, "sample")
            self.assertEqual(list(cache.iterdir()), [])
            entry.write_bytes(b">example.1\nA\n")
            with patch("urllib.request.urlopen", return_value=io.BytesIO(reference)) as fetch:
                self.assertEqual(human_reference.chromosome(cache, "sample"), b"ACGT")
                self.assertEqual(human_reference.chromosome(cache, "sample"), b"ACGT")
                fetch.assert_called_once()
            self.assertEqual(entry.read_bytes(), reference)

    def test_edits_and_reproducible_subsets(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            generate(path, [101, 150], [Decimal(".10")], 2, 23)
            with (path / "manifest.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 4)
            first_path = path / rows[0]["path"]
            before = first_path.read_bytes()
            for row in rows:
                length = int(row["length"])
                sub, ins, dele = (int(row[k]) for k in ("substitutions", "insertions", "deletions"))
                self.assertEqual(sub + ins + dele, round(length * .10))
                records = (path / row["path"]).read_text().split(">")[1:]
                ref, mut = ["".join(record.splitlines()[1:]) for record in records]
                self.assertEqual(len(ref), length)
                self.assertEqual(len(mut), length + ins - dele)
                self.assertEqual(int(row["mutant_length"]), len(mut))
                self.assertTrue(set(ref + mut) <= set("ACGT"))
            generate(path, [101], [Decimal(".1")], 1, 23)
            self.assertEqual(first_path.read_bytes(), before)


class TreeTests(unittest.TestCase):
    def test_upgma_weights_by_taxon_count(self):
        # Average linkage gives the last height (8+8+5)/6 = 3.5.
        matrix = [[0, 2, 4, 8], [2, 0, 4, 8], [4, 4, 0, 5], [8, 8, 5, 0]]
        root = build_tree(list("ABCD"), matrix, "upgma")
        self.assertEqual(root.height, 3.5)
        self.assertEqual(root.newick(), "(D:3.5,(C:2,(A:1,B:1):1):1.5)")

    def test_nj_recovers_additive_distances(self):
        matrix = [[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7],
                  [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]]
        labels = list("ABCDE")
        root = build_tree(labels, matrix, "nj")
        paths = {}

        def visit(node, path):
            if not node.children:
                paths[node.label] = path
            for child, length in node.children:
                visit(child, {**path, id(child): length})

        visit(root, {})
        for i, j in combinations(range(5), 2):
            a, b = paths[labels[i]], paths[labels[j]]
            length = sum(a[k] for k in a.keys() - b.keys()) + sum(b[k] for k in b.keys() - a.keys())
            self.assertAlmostEqual(length, matrix[i][j])

    def test_double_quotes_in_literal_tsv_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "matrix.tsv"
            path.write_text('id\t"abc"\tabc\n"abc"\t0\t1\nabc\t1\t0\n')
            labels, matrix = read_matrix(path)
            self.assertEqual(labels, ['"abc"', 'abc'])
            self.assertIn("""'"abc"'""", build_tree(labels, matrix, "nj").newick())

    def test_row_labels_and_newick_quoting(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "matrix.tsv"
            path.write_text("id\tx_y\tA'B\nA'B\t3\t0\nx_y\t0\t3\n")
            labels, matrix = read_matrix(path)
            self.assertEqual(build_tree(labels, matrix, "nj").newick(), "('x_y':1.5,'A''B':1.5)")
            path.write_text("id\ta\tb\na\t0\t2\nb\t3\t0\n")
            with self.assertRaisesRegex(ValueError, "symmetric"):
                read_matrix(path)


if __name__ == "__main__":
    unittest.main()
