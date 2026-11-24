"""Observe the unmodified window estimate and the first out-of-band traceback lookup."""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile

def replace_once(src, old, new):
    if src.count(old) != 1:
        raise ValueError("Upstream code changed")
    return src.replace(old, new)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("build", type=Path)
    p.add_argument("inputs", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--cc", default="cc")
    p.add_argument("--cxx", default="c++")
    args = p.parse_args()
    root = args.build.resolve() / "_deps/quicked-src"
    inc = [str(root), str(root / "quicked"), str(root / "quicked/include"), str(root / "quicked_utils/include")]
    q = (root / "quicked/src/quicked.c").read_text()
    q = replace_once(q, "int64_t score = windowed_matrix.cigar->score;",
                     'int64_t score = windowed_matrix.cigar->score; fprintf(stderr,"fast %lld %llu\\n",'
                     '(long long)score,(unsigned long long)windowed_matrix.high_error_window);')
    q = replace_once(q, "timer_start(aligner->timer_align);",
                     'fprintf(stderr,"final %lld\\n",(long long)score); timer_start(aligner->timer_align);')
    b = (root / "quicked/src/bpm_banded.c").read_text()
    b = replace_once(b, "int64_t h = text_length - 1;", "int reported = 0; int64_t h = text_length - 1;")
    loc = "        const uint64_t mask = 1ULL << (effective_v % BPM_W64_LENGTH);"
    b = replace_once(b, loc, """
        if (!reported && (effective_v < 0 || effective_v_r < 0 ||
                   effective_v >= 64*(int64_t)num_words64 || effective_v_r >= 64*(int64_t)num_words64)) {
            fprintf(stderr, "edge %lld %lld %lld %lld %lld\\n", (long long)h, (long long)v,
                    (long long)effective_v, (long long)effective_v_r, (long long)num_words64);
            reported = 1;
        }
""" + loc)
    with tempfile.TemporaryDirectory() as temp:
        d = Path(temp)
        for name, src in (("bound", q), ("banded", b)):
            (d / (name + ".c")).write_text(src)
            subprocess.run([args.cc, "-O3", "-march=native", *["-I"+x for x in inc],
                            "-c", d/(name+".c"), "-o", d/(name+".o")], check=True)
        subprocess.run([args.cxx, "-O3", "-march=native", "-std=c++20", "-DCHPS_BENCH_STOCK",
                        "-Iinclude", "bench/quicked/warm.cpp", "src/io/fasta.cpp", d/"bound.o", d/"banded.o",
                        args.build/"libquicked_static.a", args.build/"libedlib_static.a", "-lm", "-o", d/"run"], check=True)
        with args.output.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow("case fast_score high_error_windows final_cutoff first_h first_v relative_v relative_v_right stored_words score".split())
            for name in ("periodic512", "periodic256-tail", "identical512"):
                r = subprocess.run([d/"run", args.inputs/(name+".fa"), "quicked-stock", "1"],
                                   capture_output=True, check=True, text=True)
                diag = {fields[0]: fields[1:] for line in r.stderr.splitlines()
                        if (fields := line.split()) and fields[0] in ("fast","final","edge")}
                w.writerow([name, *diag["fast"], *diag["final"], *diag.get("edge", [""]*5),
                            r.stdout.split(",")[1]])

if __name__ == "__main__":
    main()
