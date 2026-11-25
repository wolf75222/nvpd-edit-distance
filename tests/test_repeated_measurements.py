"""Check provenance and paired blocks in the repeated and MPI benchmarks."""
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

    def test_mpi_summaries_reproduce_in_any_directory(self):
        for case in ("scaling", "collectives"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as tmp:
                archive = ROOT / "results/mpi" / case
                destination = Path(tmp)
                (destination / "trials.csv").write_bytes((archive / "trials.csv").read_bytes())
                subprocess.run([sys.executable, ROOT / "bench/mpi/summarize.py", destination], check=True)
                self.assertEqual((destination / "summary.csv").read_bytes(),
                                 (archive / "summary.csv").read_bytes())

    def test_invalid_mpi_blocks_preserve_previous_summary(self):
        archive = ROOT / "results/mpi/scaling"
        with (archive / "trials.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        mixed = [dict(row) for row in rows]
        mixed[0]["input_sha256"] = "0" * 64
        missing = [r for r in rows if (r["case"], r["block"]) != ("square", "7")]
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp)
            expected = (archive / "summary.csv").read_bytes()
            summary = destination / "summary.csv"
            summary.write_bytes(expected)
            for data, message in [(mixed, "Input or reference"), (missing, "Incomplete launch blocks")]:
                with self.subTest(message=message):
                    with (destination / "trials.csv").open("w", newline="") as source:
                        writer = csv.DictWriter(source, fieldnames=rows[0])
                        writer.writeheader()
                        writer.writerows(data)
                    run = subprocess.run([sys.executable, ROOT / "bench/mpi/summarize.py", destination],
                                         capture_output=True, text=True)
                    self.assertNotEqual(run.returncode, 0)
                    self.assertIn(message, run.stderr)
                    self.assertEqual(summary.read_bytes(), expected)
