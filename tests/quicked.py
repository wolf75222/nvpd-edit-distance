"""Compare the isolated QuickEd binary with the global edlib score."""
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main():
    quicked, editdist = sys.argv[1:]
    rng = random.Random(14)
    pairs = [("", ""), ("", "A"), ("A", "G"), ("ACGT", "ACGT"), ("ACGT", "CGTAC"),
             ("A"*150, "C"*164), ("A"*1024, "A"), ("G", "C"*480)]
    for length in (7, 31, 63, 64, 65, 127, 129, 255, 511, 2048):
        a = "".join(rng.choices("ACGT", k=length))
        pairs.append((a, "T" + a[1:-1] + "AC"))
        pairs.append((a, "".join(rng.choices("ACGT", k=length + 3))))
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "pair.fa"
        for a, b in pairs:
            path.write_text(f">a\n{a}\n>b\n{b}\n")
            ref = json.loads(subprocess.check_output([editdist, "pair", "--algorithm", "edlib", "--input", path]))
            result = json.loads(subprocess.check_output([quicked, path], timeout=15))
            if result["distance"] != ref["distance"] or result["verified_by"] != "edlib":
                raise ValueError((len(a), len(b), result, ref))
        # These periodic inputs exposed both an impossible and a plausible
        # suboptimal CIGAR score on QuickEd 1.2 with its default parameters.
        left, right = "A"*65 + "C"*64, "C"*64 + "A"*65
        for n, padding in [(512, ""), (256, "G"*256)]:
            a, b = (left*4)[:n] + padding, (right*4)[:n] + padding
            path.write_text(f">a\n{a}\n>b\n{b}\n")
            ref = json.loads(subprocess.check_output([editdist, "pair", "--algorithm", "seq", "--input", path]))
            if ref["distance"] != 128:
                raise ValueError("periodic test oracle changed")
            run = subprocess.run([quicked, path], text=True, capture_output=True, timeout=15)
            if run.returncode == 0:
                data = json.loads(run.stdout)
                if data["distance"] != 128 or data["verified_by"] != "edlib":
                    raise ValueError("accepted a suboptimal periodic alignment")
            elif "QuickEd returned " not in run.stderr or "edlib returned 128" not in run.stderr or run.stdout:
                raise ValueError(f"missing verification diagnostic: {run.stderr}")
        path.write_text(">a\nAACN\n>b\nAACG\n")
        if subprocess.run([quicked, path], capture_output=True).returncode == 0:
            raise ValueError("ambiguous input accepted")
    print(f"{len(pairs)} QuickEd scores matched edlib, 2 periodic cases checked for verified output")


if __name__ == "__main__":
    main()
