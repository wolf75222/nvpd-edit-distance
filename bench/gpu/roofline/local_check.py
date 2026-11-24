"""Verify the public GPU pair/trace and workspace matrix against edlib."""

from pathlib import Path
import json
import subprocess
import sys


def main():
    exe = Path(sys.argv[1])
    fasta = Path(sys.argv[2]) / "test.fa"
    sequences = ["", "R" * 17, "A" * 65 + "R" + "ACGT" * 70, "A" * 513 + "R", "ACGTTGCA" * 258]

    def run(mode, method, *extra):
        return subprocess.check_output(
            [exe, mode, "--algorithm", method, "--input", fasta, *extra], text=True
        )

    fasta.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(sequences)))
    if run("matrix", "nvpd-cuda-local") != run("matrix", "edlib"):
        raise ValueError("GPU local matrix differs from reference")
    for a, b in zip(sequences, sequences[1:]):
        fasta.write_text(f">a\n{a}\n>b\n{b}\n")
        opt = json.loads(run("pair", "nvpd-cuda-local", "--trace"))
        gold = json.loads(run("pair", "seq", "--trace"))
        if (opt["distance"], opt["cigar"]) != (gold["distance"], gold["cigar"]):
            raise ValueError("public traceback mismatch")
    print("10 matrix pairs and four public CUDA traces match the compiled reference")


if __name__ == "__main__":
    main()
