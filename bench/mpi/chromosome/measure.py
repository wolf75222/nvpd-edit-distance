"""Compare matrix solvers at equal core counts. The timed region includes all target indexes."""
import csv
import hashlib
import json
import random
import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mpi_benchmark import launch, run_rank_zero


def examine(matrix, lengths):
    rows = [line.split("\t") for line in matrix.splitlines()]
    labels = list(lengths)
    if rows[0] != ["id", *labels] or [r[0] for r in rows[1:]] != labels:
        raise ValueError("matrix labels or order disagree")
    values = [[int(x) for x in r[1:]] for r in rows[1:]]
    fractions = []
    eligible = certified = 0
    for i, label in enumerate(labels):
        if len(values[i]) != len(labels) or values[i][i] != 0:
            raise ValueError("invalid matrix dimensions or diagonal")
        for j in range(i):
            d = values[i][j]
            a, b = lengths[label], lengths[labels[j]]
            if d != values[j][i] or not abs(a-b) <= d <= max(a,b):
                raise ValueError("matrix violates a metric bound")
            if max(a, b) // 3 >= min(a, b):
                eligible += 1
                certified += d == abs(a-b)
            fractions.append(d/max(a,b))
    return median(fractions), eligible, certified


def configurations():
    algorithms = ["nvpd-omp", "nvpd-indexed", "edlib", "wfa2"]
    configs = [(k, nodes, a, "hilbert-weighted") for k, nodes in ((544,1),(768,1),(768,2))
               for a in algorithms]
    configs.extend((768, 2, "nvpd-indexed", policy) for policy in ("static", "steal"))
    configs = [(*c, False) for c in configs]
    configs.extend((768, 2, a, "hilbert-weighted", True) for a in ("nvpd-indexed", "edlib", "wfa2"))
    return configs


def run(binary, data, output, tmp):
    sequences = (data / "chr1.fa").read_text().splitlines()
    pairs = list(zip(sequences[::2], sequences[1::2]))
    configs = configurations()
    order = []
    for block in range(1, 6):
        random.Random(2207 + block).shuffle(configs)
        order.extend((block, *c) for c in configs)
    inputs, lengths = {}, {}
    for size in (544,768):
        count = {}
        chosen = []
        for name, seq in pairs:
            n = len(seq)
            if count.get(n,0) < size//4:
                chosen.append((name,seq))
                count[n] = count.get(n,0)+1
        inputs[size] = tmp / f"chr1-{size}.fa"
        inputs[size].write_text("".join(f"{name}\n{seq}\n" for name,seq in chosen))
        lengths[size] = {name[1:]:len(seq) for name,seq in chosen}
        if len(lengths[size]) != size: raise ValueError("incomplete set of loci")
    verified = {}
    with (output/"trials.csv").open("x") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(("records","nodes","algorithm","schedule","length_bound","trial","cores","pairs","equivalent_cells",
                    "seconds","max_rank_rss_kib","stolen","cross_node","certified","median_relative_distance","matrix_sha256"))
        for block,size,nodes,algorithm,policy,certificate in order:
            args = [binary, "matrix", "--input", inputs[size], "--schedule", policy, "--algorithm", algorithm,
                    "--threads", "8", "--grain", "4", "--warmups","1","--trials","1","--report-rss","1"]
            if certificate:
                args.extend(["--length-bound", "1"])
            proc = run_rank_zero(launch(nodes,24,8), args, tmp)
            digest = hashlib.sha256(proc.stdout.encode()).hexdigest()
            if size in verified:
                expected,relative,oracle_count = verified[size]
                if expected != digest: raise ValueError(f"{algorithm} disagrees on {size} loci")
            else:
                relative, _, oracle_count = examine(proc.stdout, lengths[size])
                verified[size] = digest,relative,oracle_count
            diags = [json.loads(line) for line in proc.stderr.splitlines() if line.startswith("{")]
            timing = [r for r in diags if "mode" in r]
            if len(timing) != 1: raise ValueError("incomplete trial series")
            seq_lengths = list(lengths[size].values())
            area = (sum(seq_lengths)**2-sum(x*x for x in seq_lengths))//2
            for t in timing:
                ranks = [r for r in diags if r.get("trial") == t["trial"] and "rank" in r]
                npairs = sum(r["pairs"] for r in ranks)
                reported = sum(r["certified_pairs"] for r in ranks)
                if (sorted(r["rank"] for r in ranks) != list(range(nodes*24))
                        or t["ranks"] != nodes*24 or t["threads_per_rank"] != 8):
                    raise ValueError("missing rank or worker-team diagnostics")
                if reported != (oracle_count if certificate else 0):
                    raise ValueError(f"{algorithm} reported {reported} certificates; expected "
                                     f"{oracle_count if certificate else 0}")
                if (t["algorithm"] != algorithm or t["sampled_physical_cores"] != 192*nodes
                    or npairs != size*(size-1)//2 or sum(r["cells"] for r in ranks) != area):
                    raise ValueError(f"incomplete allocation or matrix: {t}, {npairs}")
                w.writerow((size,nodes,algorithm,policy,int(certificate),block,t["sampled_physical_cores"],npairs,area,
                            t["seconds"],t["max_rank_rss_kib"],sum(r["stolen_pairs"] for r in ranks),
                            sum(r["cross_node_pairs"] for r in ranks),reported,relative,digest))
            out.flush()
            print(size,nodes,algorithm,policy,"bound" if certificate else "direct",timing[-1]["seconds"],flush=True)


if __name__ == "__main__":
    run(*(Path(arg).resolve() for arg in sys.argv[1:]))
