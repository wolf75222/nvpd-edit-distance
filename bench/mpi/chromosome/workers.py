"""Compare indexed NvPD on the same EPYC cores under OpenMP and Kokkos OpenMP."""
import csv
import hashlib
import json
import random
import statistics
import sys
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mpi_benchmark import launch, run_rank_zero

POLICIES = ("hilbert-weighted", "steal")
BACKENDS = ("openmp", "kokkos")
BLOCKS = range(1, 6)


def summarize(out, rows):
    keyed = {(b, p, n, a): t for b, p, n, a, t in rows}
    if set(keyed) != set(product(BLOCKS, POLICIES, (1, 2), BACKENDS)) or len(rows) != len(keyed):
        raise ValueError("incomplete or duplicate matrix configurations")
    with (out / "summary.csv").open("x") as stream:
        wr = csv.writer(stream, lineterminator="\n")
        wr.writerow(("schedule", "nodes", "backend", "metric", "blocks", "mean", "sd", "low95", "high95"))
        for p, n, a in product(POLICIES, (1, 2), BACKENDS):
            times = [keyed[b,p,n,a] for b in BLOCKS]
            metrics = {"seconds": times}
            if n == 2:
                speedup = [keyed[b,p,1,a] / keyed[b,p,2,a] for b in BLOCKS]
                metrics.update(speedup=speedup, efficiency=[s / 2 for s in speedup])
            if a == "kokkos":
                metrics["openmp_over_kokkos"] = [keyed[b,p,n,"openmp"] / keyed[b,p,n,a] for b in BLOCKS]
            for metric, values in metrics.items():
                mean, sd = statistics.mean(values), statistics.stdev(values)
                ci = 2.7764451052 * sd / len(values)**.5
                wr.writerow((p,n,a,metric,len(values),mean,sd,mean-ci,mean+ci))


def run(binary, fasta, reference, out, tmp):
    oracle = json.loads(reference.read_text())
    if hashlib.sha256(fasta.read_bytes()).hexdigest() != oracle["fasta_sha256"]:
        raise ValueError("FASTA differs from the edlib reference")
    sequences = fasta.read_bytes().splitlines()
    sizes = [len(s) for s in sequences[1::2]]
    if len(sizes) != oracle["records"] or any(not s.startswith(b">") for s in sequences[::2]):
        raise ValueError("expected the two-line FASTA of the chromosome cohort")
    cells = (sum(sizes)**2 - sum(s*s for s in sizes)) // 2
    configs = list(product(POLICIES, (1, 2), BACKENDS))
    samples = []
    with (out / "trials.csv").open("x") as stream:
        wr = csv.writer(stream, lineterminator="\n")
        wr.writerow(("block", "schedule", "nodes", "backend", "physical_cores", "active_workers",
                     "seconds", "max_rank_rss_kib", "pairs", "equivalent_cells", "stolen_pairs",
                     "cross_node_pairs", "matrix_sha256"))
        for b in BLOCKS:
            order = configs.copy()
            random.Random(4319+b).shuffle(order)
            for policy, nodes, backend in order:
                args = [binary, "matrix", "--input", fasta, "--algorithm", "nvpd-indexed", "--schedule", policy,
                        "--workers", backend, "--threads", "8", "--grain", "4", "--warmups", "1",
                        "--trials", "1", "--report-rss", "1"]
                proc = run_rank_zero(launch(nodes,24,8), args, tmp)
                digest = hashlib.sha256(proc.stdout.encode()).hexdigest()
                diagnostics = [json.loads(line) for line in proc.stderr.splitlines() if line.startswith("{")]
                timing, = (r for r in diagnostics if "mode" in r)
                ranks = [r for r in diagnostics if "rank" in r]
                npairs = sum(r["pairs"] for r in ranks)
                ncores = 192 * nodes
                if (digest != oracle["matrix_sha256"] or npairs != oracle["pairs"]
                    or sorted(r["rank"] for r in ranks) != list(range(24*nodes))
                    or timing["worker_backend"] != backend or timing["algorithm"] != "nvpd-indexed"
                    or timing["ranks"] != 24*nodes or timing["threads_per_rank"] != 8
                    or timing["sampled_physical_cores"] != ncores or not 0 < timing["assigned_workers"] <= ncores
                    or sum(r["cells"] for r in ranks) != cells or any(r["certified_pairs"] for r in ranks)):
                    raise ValueError(f"score, placement or pair assignment differs: {timing}")
                samples.append((b, policy, nodes, backend, timing["seconds"]))
                wr.writerow((b,policy,nodes,backend,ncores,timing["assigned_workers"],timing["seconds"],
                             timing["max_rank_rss_kib"],npairs,cells,sum(r["stolen_pairs"] for r in ranks),
                             sum(r["cross_node_pairs"] for r in ranks),digest))
                stream.flush()
    summarize(out, samples)


if __name__ == "__main__":
    run(*(Path(s).resolve() for s in sys.argv[1:]))
