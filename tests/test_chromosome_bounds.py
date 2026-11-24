import gzip
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
from chromosome_bounds import bounds
from fasta_records import read_pair


def scalar(a, b):
    prev = list(range(len(b)+1))
    for i, x in enumerate(a, 1):
        row = [i]
        for j, y in enumerate(b, 1):
            row.append(min(prev[j]+1, row[-1]+1, prev[j-1]+(x != y)))
        prev = row
    return prev[-1]


class ChromosomeBounds(unittest.TestCase):
    def test_constructive_bounds(self):
        rng = random.Random(23)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "chain.gz")
            with gzip.open(path, "wt") as out:
                out.write("chain 1 chr21 16 + 1 14 chr21 19 + 2 17 23\n5 2 4\n6\n")
            for _ in range(80):
                a, b = [bytes(rng.choices(b"ACGTN", k=n)) for n in (16,19)]
                r = bounds(a, b, path, 23)
                distance = scalar(a, b)
                self.assertLessEqual(r["lower"], distance)
                self.assertLessEqual(distance, r["upper"])

    def test_malformed_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "chain.gz")
            for header, blocks in [("0 4", "4 0 0"), ("0 5", "5"), ("0 4", "4 -1 -1\n1")]:
                with gzip.open(path, "wt") as out:
                    out.write(f"chain 1 chr21 4 + {header} chr21 4 + {header} 23\n{blocks}\n")
                with self.assertRaises(ValueError):
                    bounds(b"ACGT", b"ACGT", path, 23)

    def test_pair_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "input.fa")
            path.write_bytes(b"\n>a desc > example\na c\tG\n>b\nnr\n")
            self.assertEqual(read_pair(path), [b"ACG", b"NR"])
            for invalid in [b"ACT\n>b\nA\n", b">a\nA\n>a\nG\n", b">a\nA?\n>b\nG\n"]:
                path.write_bytes(invalid)
                with self.assertRaises(ValueError):
                    read_pair(path)
