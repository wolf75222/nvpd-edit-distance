"""Exercise FASTA normalization and malformed-input errors through the CPU CLI."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "input.fa"

        def read(contents):
            path.write_bytes(contents)
            return subprocess.run([binary, "pair", "--algorithm", "edlib", "--input", path],
                                  capture_output=True, text=True)

        for width in (0, 1, 15, 16, 31, 32, 63, 64, 65, 10_001):
            canonical = (b"ACGT" * ((width + 3) // 4))[:width]
            letters = b"RYSWKMBDHVN"
            sample = (b" \t\r\n>first description\r\n" + canonical + b"\n" +
                      b" \t".join(bytes((x,)) for x in letters.lower()) + b"\v\f\n" +
                      b">second\n" + canonical.lower() + letters + b"\n")
            result = read(sample)
            if result.returncode:
                raise AssertionError(result.stderr)
            row = json.loads(result.stdout)
            if (row["distance"], row["m"], row["n"]) != (0, width + len(letters), width + len(letters)):
                raise AssertionError(row)
        for malformed, message in ((b">a\nAC\x00T\n>b\nA", "invalid DNA symbol"),
                                   (b">a\nA\xffT\n>b\nA", "invalid DNA symbol"),
                                   (b" \tA\n>a\n", "sequence before FASTA header"),
                                   (b">a desc\nA\n>a other\nA", "duplicate FASTA ID")):
            result = read(malformed)
            if result.returncode == 0 or message not in result.stderr:
                raise AssertionError((result.returncode, result.stderr))
    print("10 line widths and 4 malformed FASTA cases passed")


if __name__ == "__main__":
    main()
