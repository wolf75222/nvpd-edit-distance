"""Read a complete-genome package from the NCBI Datasets CLI."""
import json
import zipfile


def read_package(package, tax_id):
    with zipfile.ZipFile(package) as archive:
        base = "ncbi_dataset/data/"
        metadata = {}
        for line in archive.open(base + "data_report.jsonl"):
            row = json.loads(line)
            if row["virus"]["taxId"] == tax_id and row["completeness"] == "COMPLETE":
                if row["accession"] in metadata:
                    raise ValueError(f"duplicate NCBI accession {row['accession']}")
                metadata[row["accession"]] = row
        sequences = {}
        with archive.open(base + "genomic.fna") as fasta:
            for raw in fasta:
                line = raw.decode("ascii").strip()
                if line.startswith(">"):
                    accession = line[1:].split()[0]
                    if accession in sequences:
                        raise ValueError(f"duplicate FASTA record {accession}")
                    sequences[accession] = []
                elif line:
                    sequences[accession].append(line)
    return metadata, {key: "".join(parts).upper() for key, parts in sequences.items()}
