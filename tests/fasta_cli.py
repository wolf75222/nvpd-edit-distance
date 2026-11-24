"""Exercise FASTA normalization and malformed-input errors through the CPU CLI."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def vector_boundaries(binary, path, index):
    canonical = b"ACGT" * 16
    # Bytes 15/16 fall on opposite sides of the SIMD load boundary.
    for pos, symbol in zip((0, 15, 16, 31), b"rybn"):
        changed = canonical[:pos] + bytes([symbol]) + canonical[pos + 1:]
        for mode in ("pair", "matrix"):
            records = (changed, changed.upper()) if mode == "pair" else (changed, changed.upper(), canonical)
            raw = bytearray()
            entries = []
            for name, bases in zip(("a", "b", "c"), records):
                raw.extend(f">{name}\n".encode())
                entries.append(f"{name}\t64\t{len(raw)}\t64\t65\n")
                raw.extend(bases + b"\n")
            index.write_text("".join(entries))
            for args in ([], ["--fai", index]):
                argv = [binary, mode, "--algorithm", "edlib", "--input", path, *args]
                path.write_bytes(raw)
                output = subprocess.check_output(argv, text=True)
                if mode == "pair":
                    row = json.loads(output)
                    if (row["m"], row["n"], row["distance"]) != (64, 64, 0):
                        raise AssertionError((pos, row))
                elif output != "id\ta\tb\tc\na\t0\t0\t1\nb\t0\t0\t1\nc\t1\t1\t0\n":
                    raise AssertionError((pos, output))
                for invalid in (0, 255):
                    bad = raw.copy()
                    bad[3 + pos] = invalid  # first base follows ">a\n"
                    path.write_bytes(bad)
                    result = subprocess.run(argv, capture_output=True, text=True)
                    if not result.returncode or "invalid DNA symbol in a" not in result.stderr:
                        raise AssertionError((pos, invalid, result.stderr))


def main():
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "input.fa"

        def read(contents, *args):
            path.write_bytes(contents)
            return subprocess.run([binary, "pair", "--algorithm", "edlib", "--input", path, *args],
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
        index = path.with_suffix(".fa.fai")
        # Each record can have its own width and line terminator.
        sample = b">a description\r\nacgt\r\nRyN\r\n>b\nACG\nTRY\nN"
        a, b = sample.index(b"acgt"), sample.index(b"ACG")
        row_a = f"a\t7\t{a}\t4\t6\n"
        row_b = f"b\t7\t{b}\t3\t4\n"
        index.write_text(row_a + row_b)
        indexed = read(sample, "--fai", index)
        if indexed.returncode or json.loads(indexed.stdout)["distance"] != 0:
            raise AssertionError(indexed.stderr or indexed.stdout)
        matrix = subprocess.run([binary, "matrix", "--algorithm", "edlib", "--input", path, "--fai", index],
                                check=True, capture_output=True, text=True)
        plain = subprocess.check_output([binary, "matrix", "--algorithm", "edlib", "--input", path], text=True)
        if matrix.stdout != plain:
            raise AssertionError("indexed matrix differs")
        invalid = (row_a, row_b + row_a, row_a + row_b + row_b, row_a.replace("a", "wrong") + row_b,
                   row_a.replace(f"\t{a}\t", f"\t{a+1}\t") + row_b, row_a.replace("7", "6") + row_b,
                   row_a.replace("7", "8") + row_b, row_a.replace("4\t6", "4\t5") + row_b,
                   row_a.replace("7", "-7") + row_b, row_a.replace("7", "18446744073709551616") + row_b,
                   row_a.rstrip() + "\textra\n" + row_b)
        for contents in invalid:
            index.write_text(contents)
            result = read(sample, "--fai", index)
            if result.returncode == 0 or "FAI" not in result.stderr:
                raise AssertionError((contents, result.stderr))
        index.unlink()
        result = read(sample, "--fai", index)
        if result.returncode == 0 or "cannot open FAI" not in result.stderr:
            raise AssertionError(result.stderr)
        empty = b">a\n>b\nA"
        index.write_text("a\t0\t3\t0\t0\nb\t1\t6\t1\t1\n")
        result = read(empty, "--fai", index)
        if result.returncode or json.loads(result.stdout)["distance"] != 1:
            raise AssertionError(result.stderr or result.stdout)
        # Samtools 1.21 omits the leading empty record and reports LINEWIDTH=2
        # for this one-base final line even when the FASTA lacks its final LF.
        index.write_text("b\t1\t6\t1\t2\n")
        result = read(empty, "--fai", index)
        if result.returncode or json.loads(result.stdout)["distance"] != 1:
            raise AssertionError(result.stderr or result.stdout)
        # Match across several omitted records and keep rejecting any omitted
        # nonempty record, including after the last FAI line.
        cohort = b">a\nA\n>b\n>c\n>d\nCC\n"
        index.write_text(f"a\t1\t3\t1\t2\nd\t2\t{cohort.index(b'CC')}\t2\t3\n")
        path.write_bytes(cohort)
        matrix_args = [binary, "matrix", "--algorithm", "edlib", "--input", path]
        if subprocess.check_output(matrix_args + ["--fai", index]) != subprocess.check_output(matrix_args):
            raise AssertionError("matrix disagrees across omitted empty FASTA records")
        for invalid in (cohort.replace(b">b\n", b">b\nA\n"), cohort + b">e\nA\n"):
            path.write_bytes(invalid)
            result = subprocess.run(matrix_args + ["--fai", index], capture_output=True, text=True)
            if not result.returncode or "FAI" not in result.stderr:
                raise AssertionError("indexed reader accepted missing nonempty record")
        vector_boundaries(binary, path, index)
    print("FASTA normalization, indexed pair/matrix, empty records, and stale/malformed index checks passed")


if __name__ == "__main__":
    main()
