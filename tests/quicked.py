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
    pairs = [("", ""), ("", "A"), ("A", ""), ("A", "G"), ("ACGT", "ACGT"), ("ACGT", "CGTAC"),
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
            verifier = "edlib" if a and b else "length"
            if result["distance"] != ref["distance"] or result["verified_by"] != verifier:
                raise ValueError((len(a), len(b), result, ref))
            if verifier == "length" and (result["source"], result["retries"]) != ("empty-input", 0):
                raise ValueError(result)
        # These periodic inputs exposed both an impossible and a plausible
        # suboptimal CIGAR score on QuickEd 1.2 with its default parameters.
        left, right = "A"*65 + "C"*64, "C"*64 + "A"*65
        witnesses = [((left*4)[:n] + tail, (right*4)[:n] + tail, 128)
                     for n, tail in [(512, ""), (256, "G"*256)]]
        x = ("A"*64 + "C"*64)*2
        witnesses.append((x, x.translate(str.maketrans("AC", "CA")), 128))
        x = "A"*64 + "C"*128 + "A"*128 + "C"*128 + "A"*64 + "G"*256
        witnesses.append((x, x.translate(str.maketrans("AC", "CA")), 256))
        for a, b, expected in witnesses:
            path.write_text(f">a\n{a}\n>b\n{b}\n")
            ref = json.loads(subprocess.check_output([editdist, "pair", "--algorithm", "seq", "--input", path]))
            if ref["distance"] != expected:
                raise ValueError("periodic test oracle changed")
            data = json.loads(subprocess.check_output([quicked, path], timeout=15))
            if (data["distance"], data["first_score"], data["source"], data["retries"]) != (
                    expected, expected, "quicked-default", 0):
                raise ValueError(f"boundary guard failed on periodic alignment: {data}")
        with path.open("rb") as unwritable:
            result = subprocess.run([quicked, path], stdout=unwritable, stderr=subprocess.PIPE)
        if not result.returncode or b"cannot write result" not in result.stderr:
            raise ValueError("QuickEd did not report its output error")
        path.write_text(">a\nAACN\n>b\nAACG\n")
        if subprocess.run([quicked, path], capture_output=True).returncode == 0:
            raise ValueError("ambiguous input accepted")
    print(f"{len(pairs)} QuickEd scores matched edlib, 4 periodic regression cases passed")


if __name__ == "__main__":
    main()
