"""Check that matched benchmark blocks use the same executable."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class RepeatedMeasurements(unittest.TestCase):
    def test_mixed_executable_rejected_before_summary(self):
        with (ROOT / "results/repeated/epyc-trace/trials.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            fields = reader.fieldnames
            rows = list(reader)
        rows[-1]["binary_sha256"] = "f" * 64
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp)
            with (raw / "lane-0.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fields)
                writer.writeheader()
                writer.writerows(rows)
            result = subprocess.run(
                [sys.executable, ROOT / "bench/repeated/collect.py", "epyc-trace", raw, raw / "out"],
                capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("inconsistent executable hash", result.stderr)
            self.assertFalse((raw / "out").exists())
