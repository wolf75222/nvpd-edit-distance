import csv
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

STUDY = Path(__file__).resolve().parents[1] / "bench/human/grace-numa/amortization"
with patch.object(sys, "path", [str(STUDY), *sys.path]):
    import analyze_e2e
    import measure_e2e


class GraceAmortization(unittest.TestCase):
    def test_design_balance_and_target_grid(self):
        for n in (12, 16):
            labels = measure_e2e.balanced(n, 7)
            self.assertEqual(len(labels), n)
            for r in range(4):
                self.assertEqual(sum(rotation == r for rotation, _ in labels), n // 4)
            self.assertEqual(sum(reverse for _, reverse in labels), n // 2)
            if n == 16:
                for label in labels:
                    self.assertEqual(labels.count(label), 2)
        with self.assertRaisesRegex(ValueError, "Invalid completed"):
            measure_e2e.timing_row("timing 7 1 1 0 1 1 6".split(), 7)

    def test_intervals_include_zero_and_preserve_parity(self):
        self.assertAlmostEqual(analyze_e2e.student_quantile(15, .975), 2.13144954556, places=8)
        self.assertLess(analyze_e2e.estimates([1., -1.] * 8, 10)["family_low95"],
                        analyze_e2e.estimates([1., -1.] * 8, 10)["low95"])
        data = {k: dict(mean=v, family_low95=v-.25, family_high95=v+.25)
                for k, v in [(1, -2.), (2, .5), (3, -.1), (4, 1.)]}
        result = analyze_e2e.thresholds(data, "direct", "gain_s", "all", sorted(data))
        self.assertEqual((result["first_mean_k"], result["sustained_mean_k"],
                          result["last_family_loss_k"], result["sustained_family_gain_k"]), (2, 4, 1, 4))
        self.assertEqual(analyze_e2e.thresholds(data, "direct", "gain_s", "odd", sorted(data))["sustained_mean_k"], "")

    def test_paired_output_and_incomplete_confirmation(self):
        rows = []
        def add(phase, k, b, p, assignment):
            rotation, reverse = assignment
            for pos, m in enumerate(measure_e2e.MODES[::-1] if reverse else measure_e2e.MODES, 1):
                read = 1.
                copy = .075 if m == "local" else 0.
                align = k * (.01 if m == "local" else .02)
                total = read + copy + align
                rows.append(dict(phase=phase, block=b, order=p, position=pos, mode=m, rotation=rotation,
                                 k=k, total_s=total, read_s=read, copy_s=copy, setup_s=read+copy,
                                 align_s=align, scores=k, process_s=(total+.003 if phase == "direct" else "")))
        n = measure_e2e.PREFIX_BLOCKS
        for b, assignment in enumerate(measure_e2e.balanced(n, 0), 1):
            for k in range(2, measure_e2e.MAX_K+1, 2):
                add("prefix", k, b, 1, assignment)
        cross, keys, blocks, estimate = measure_e2e.confirm_plan(rows)
        self.assertEqual(cross, 8)
        self.assertTrue(set(range(1, 15)).issubset(keys))
        for b in range(1, blocks+1):
            order = keys[:]
            random.Random(measure_e2e.SEED + 1000 + b).shuffle(order)
            for i, k in enumerate(order, 1):
                add("direct", k, b, i, measure_e2e.balanced(blocks, k)[b-1])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            analyze_e2e.write_csv(path/"timings.csv", rows)
            analyze_e2e.write_csv(path/"confirm-design.csv", [
                dict(pilot_sustained_k=cross, direct_blocks=blocks, k=k, estimated_campaign_s=estimate) for k in keys])
            analyze_e2e.analyze(path)
            with (path/"threshold.csv").open() as f:
                target = [r for r in csv.DictReader(f) if r["phase"] == "direct" and r["parity"] == "all"]
            self.assertEqual([int(r["sustained_family_gain_k"]) for r in target], [8, 8])
            analyze_e2e.write_csv(path/"timings.csv", rows[:-1])
            with self.assertRaisesRegex(ValueError, "Missing an expected timing"):
                analyze_e2e.analyze(path)


if __name__ == "__main__":
    unittest.main()
