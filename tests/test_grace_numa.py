import csv
from collections import Counter
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

BENCH = Path(__file__).resolve().parents[1] / "bench/human/grace-numa"


def load(name):
    spec = importlib.util.spec_from_file_location(name, BENCH / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


measure = load("measure")
with patch.dict(sys.modules, {"measure": measure}):
    analyze = load("analyze")


class GraceNUMA(unittest.TestCase):
    def test_balanced_blocks(self):
        orders = measure.design()
        self.assertEqual(len(orders), 14)
        self.assertTrue(all(set(order) == set(measure.MODES) for order in orders))
        precedences = Counter(pair for order in orders for pair in zip(order, order[1:]))
        positions = Counter(pair for order in orders for pair in enumerate(order))
        self.assertEqual((len(precedences), set(precedences.values())), (42, {2}))
        self.assertEqual((len(positions), set(positions.values())), (49, {2}))

    def test_residency_maps_are_exclusive_and_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = Path(tmp)
            (proc / "maps").write_text("1000-2000 rw-p 0000 00:00 0\n2000-3000 rw-p 0000 00:00 0\n")
            (proc / "numa_maps").write_text("1000 default N0=1 kernelpagesize_kB=4\n"
                                            "2000 default N1=1 kernelpagesize_kB=4\n")
            buffers = [(0, 0, 0, 0, 4096, 5000)]
            self.assertEqual(list(measure.snapshot(proc, buffers, 4096))[0][-1], {0: 4096, 1: 4096})
            with self.assertRaisesRegex(RuntimeError, "dedicated"):
                list(measure.snapshot(proc, [(0, 0, 0, 0, 4224, 5000)], 4096))
            (proc / "numa_maps").write_text("1000 default N0=1 kernelpagesize_kB=4\n")
            with self.assertRaisesRegex(RuntimeError, "resident bytes"):
                list(measure.snapshot(proc, buffers, 4096))

    def test_contrasts_and_control_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "inputs.tsv"
            manifest.write_text("chr.fa\t5000\t5000\t0\n" * 4)
            (root / "topology.txt").write_text("base_page_bytes=4096\ncpu=0 core=0 package=0 node=0 distances=10 20\n"
                                               "cpu=72 core=1 package=1 node=1 distances=20 10\n")
            reference = dict(single=1.2, shared0=1.4, shared1=1.5, local=0.8, remote=1.6, twins0=1.3, twins1=1.4)
            with (root / "trials.csv").open("w", newline="") as times, (root / "residency.csv").open("w", newline="") as pages:
                tc, pc = csv.writer(times), csv.writer(pages)
                tc.writerow(measure.TRIALS)
                pc.writerow(measure.RESIDENCY)
                for b, order in enumerate(measure.design(), 1):
                    for i, mode in enumerate(order, 1):
                        one = mode == "single"
                        tc.writerow((b, i, mode, 1 if one else 2, 64, 640000, .05,
                                     reference[mode] + b/1000, "0" if one else "0|72"))
                        for phase in ("before", "after"):
                            for c, owner in enumerate(measure.panel_owners(mode)):
                                node = int(owner == 72)
                                for k in range(4):
                                    for side in (0, 1):
                                        pc.writerow((b, mode, phase, c, k, side, owner, node, node, 5000, 8192, 8192))
            analyze.analyze(root, manifest)
            summary = {r["metric"]: r for r in analyze.read_csv(root / "summary.csv")}
            self.assertAlmostEqual(float(summary["remote_minus_local_s"]["mean"]), .8)
            self.assertAlmostEqual(float(summary["shared0_minus_twins0_s"]["mean"]), .1)
            self.assertAlmostEqual(float(summary["twins0_minus_local_s"]["mean"]), .5)
            self.assertAlmostEqual(float(analyze.read_csv(root / "attribution.csv")[-1]["ratio"]), 3.)
            # For one panel, half of one chromosome sequence is now remote.
            residency = analyze.read_csv(root / "residency.csv")
            obs = residency[0]
            obs["resident_bytes"] = 4096
            residency.append(dict(obs, observed_node=1-int(obs["expected_node"])))
            analyze.write_csv(root / "residency.csv", residency)
            with self.assertRaisesRegex(ValueError, ">99%"):
                analyze.analyze(root, manifest)


if __name__ == "__main__":
    unittest.main()
