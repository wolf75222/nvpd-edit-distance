"""Reconstruct phased GRCh37 chromosome 21 windows from the public IGSR phase-three VCF."""
import argparse
import collections
import hashlib
import json
import random
from pathlib import Path
from urllib.request import urlopen

import pysam

VCF = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
       "ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
START, END = 20_000_000, 20_032_768
SEED = 1001


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--samples", type=int, default=512)
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    reference_url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr21;start={START};end={END}"
    with urlopen(reference_url, timeout=45) as response:
        ref = json.load(response)["dna"].upper()
    if len(ref) != END - START or set(ref) - set("ACGT"):
        raise ValueError("ambiguous or incomplete reference interval")
    index = out / "chr21.tbi"
    with urlopen(VCF + ".tbi", timeout=45) as response:
        index.write_bytes(response.read())
    vcf = pysam.TabixFile(VCF, index=str(index))
    header = next(line for line in vcf.header if line.startswith("#CHROM"))
    ids = header.split("\t")[9:]
    selected = random.Random(SEED).sample(range(len(ids)), args.samples)
    chunks = [[] for _ in range(2 * len(selected))]
    cursor = START
    counts = collections.Counter()
    variant_rows = hashlib.sha256()
    for line in vcf.fetch("21", START, END):
        fields = line.split("\t")
        pos = int(fields[1]) - 1
        src, alt = fields[3:5]
        counts["total_records"] += 1
        if pos < cursor or pos + len(src) > END:
            counts["overlap_or_boundary"] += 1
            continue
        if fields[6] != "PASS" or any(c not in "ACGT" for c in src + alt):
            counts["filtered_or_nonbiallelic"] += 1
            continue
        if src != ref[pos-START:pos-START+len(src)]:
            raise ValueError(f"reference mismatch at {pos + 1}")
        gt = fields[8].split(":").index("GT")
        calls = [fields[9+i].split(":")[gt] for i in selected]
        if any(call not in {"0|0", "1|0", "0|1", "1|1"} for call in calls):
            counts["unphased_or_missing"] += 1
            continue
        flank = ref[cursor-START:pos-START]
        for i, call in enumerate(calls):
            for h, allele in enumerate((call[0], call[2])):
                chunks[2*i+h].extend((flank, src if allele == "0" else alt))
        counts["snps" if len(src) == len(alt) == 1 else "other_small_variants"] += 1
        variant_rows.update((line + "\n").encode())
        cursor = pos + len(src)
    vcf.close()
    tail = ref[cursor-START:]
    records = ["".join(parts) + tail for parts in chunks]
    fasta = out / "cohort.fa"
    with fasta.open("w") as output:
        for i, seq in enumerate(records):
            output.write(f">hap{i:04d}\n{seq}\n")
    report = dict(reference_url=reference_url, reference_sha256=hashlib.sha256(ref.encode()).hexdigest(),
                  vcf=VCF, vcf_index_sha256=digest(index), retained_rows_sha256=variant_rows.hexdigest(),
                  seed=SEED, sampled_individuals=len(selected), haplotypes=len(chunks),
                  distinct_sequences=len(set(records)), available_individuals=len(ids), counts=dict(counts),
                  length_min=min(map(len, records)), length_max=max(map(len, records)), fasta_sha256=digest(fasta))
    (out / "dataset.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
