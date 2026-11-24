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
            if result["distance"] != ref["distance"]:
                raise ValueError((len(a), len(b), result, ref))
        path.write_text(">a\nAACN\n>b\nAACG\n")
        if subprocess.run([quicked, path], capture_output=True).returncode == 0:
            raise ValueError("ambiguous input accepted")
    print(f"{len(pairs)} QuickEd scores matched edlib")


if __name__ == "__main__":
    main()
