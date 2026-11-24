"""Read paired FASTA with the executable's literal IUPAC and normalization rules."""
from pathlib import Path


def read_pair(path: Path):
    ids, sequences = set(), []
    with path.open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                header = line[1:].split()
                if not header or header[0] in ids:
                    raise ValueError(f"missing or duplicate FASTA ID in {path}")
                if len(sequences) == 2:
                    raise ValueError(f"expected two FASTA records in {path}")
                ids.add(header[0])
                sequences.append(bytearray())
            else:
                bases = b"".join(line.split()).upper()
                if not bases:
                    continue
                if not sequences or bases.translate(None, b"ACGTRYSWKMBDHVN"):
                    raise ValueError(f"invalid FASTA sequence in {path}")
                sequences[-1].extend(bases)
    if len(sequences) != 2:
        raise ValueError(f"expected two FASTA records in {path}")
    return sequences
