"""Check adjacent-pair distribution and its distinct RSS reporting option."""
import argparse
import csv
import io
import json
from pathlib import Path
import shlex
import subprocess
import tempfile


def check(binary, launcher):
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "pairs.fa"
        cmd = [*shlex.split(launcher), binary, "pairs", "--input", str(fasta), "--threads", "2"]
        fasta.write_text(">a\nACG\n>b\nATG\n>c\nACTTTT\n>d\nA\n>e\nGGGGG\n>f\nGGGGG\n")
        result = subprocess.run([*cmd,"--warmups","1","--report-rss","1"], check=True, capture_output=True, text=True)
        distances = list(csv.reader(io.StringIO(result.stdout),delimiter="\t"))
        if distances != [["id_a","id_b","distance"],["a","b","1"],["c","d","5"],["e","f","0"]]:
            raise ValueError("incorrect adjacent-pair scores or order")
        [timing] = [json.loads(line) for line in result.stderr.splitlines() if line.startswith("{")]
        if timing["max_rank_rss_kib"] <= 0 or timing["assigned_workers"] != 3:
            raise ValueError("invalid observed RSS or incomplete worker participation")
        fasta.write_text(">a\nACG\n>b\nATG\n>c\nACTTTT\n")
        if not subprocess.run(cmd,capture_output=True).returncode:
            raise ValueError("odd input record count should fail")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary")
    parser.add_argument("--launcher",default="mpiexec -n 2")
    args = parser.parse_args()
    check(args.binary,args.launcher)
