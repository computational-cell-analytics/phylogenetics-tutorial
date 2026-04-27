#!/usr/bin/env python3

"""
Download a small primate mitochondrial dataset for a phylogenetics exercise.

Dataset:
    - Human
    - Chimpanzee
    - Gorilla
    - Orangutan
    - Gibbon
    - Rhesus macaque

The script downloads complete mitochondrial RefSeq records from NCBI,
extracts the cytochrome b gene (MT-CYB / CYTB), writes a FASTA file,
computes a simple p-distance matrix, and optionally plots an UPGMA-like
average-linkage tree.

Dependencies:
    conda install -c conda-forge biopython scipy matplotlib pandas
"""

from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


# IMPORTANT: NCBI asks users of Entrez to provide an email address.
Entrez.email = "constantin.pape@inormatik.uni-goettingen.de"


ACCESSIONS = {
    "Human": "NC_012920.1",
    "Chimpanzee": "NC_001643.1",
    "Gorilla": "NC_001645.1",
    "Orangutan": "NC_002083.1",      # Pongo abelii
    "Gibbon": "NC_002082.1",         # Hylobates lar
    "Macaque": "NC_005943.1",        # Macaca mulatta
}


OUT_FASTA = Path("primate_mt_cytb.fasta")


def download_genbank_records(accessions):
    """Download GenBank records from NCBI."""
    ids = ",".join(accessions.values())
    with Entrez.efetch(
        db="nuccore",
        id=ids,
        rettype="gb",
        retmode="text",
    ) as handle:
        records = list(SeqIO.parse(handle, "genbank"))

    by_accession = {record.id: record for record in records}
    return by_accession


def feature_matches_cytb(feature):
    """Return True if a GenBank feature appears to describe cytochrome b."""
    qualifiers = feature.qualifiers

    texts = []
    for key in ("gene", "product", "note"):
        texts.extend(qualifiers.get(key, []))

    text = " ".join(texts).lower()

    return (
        "cytb" in text
        or "mt-cyb" in text
        or "cytochrome b" in text
        or "cytochrome-b" in text
    )


def extract_cytb(record, label):
    """Extract the cytochrome b sequence from a mitochondrial GenBank record."""
    candidates = []

    for feature in record.features:
        if feature.type in {"gene", "CDS"} and feature_matches_cytb(feature):
            seq = feature.location.extract(record.seq)
            candidates.append((feature.type, seq))

    if not candidates:
        raise RuntimeError(f"No CYTB / cytochrome b feature found for {label}.")

    # Prefer CDS over gene if both are present.
    candidates = sorted(candidates, key=lambda x: 0 if x[0] == "CDS" else 1)
    seq = candidates[0][1]

    return SeqRecord(
        seq,
        id=label,
        name=label,
        description=f"{label} mitochondrial cytochrome b extracted from {record.id}",
    )


def write_cytb_fasta(accessions, out_path):
    records_by_acc = download_genbank_records(accessions)

    cytb_records = []
    for label, accession in accessions.items():
        if accession not in records_by_acc:
            raise RuntimeError(f"Could not find downloaded record for {label}: {accession}")

        cytb = extract_cytb(records_by_acc[accession], label)
        cytb_records.append(cytb)

    SeqIO.write(cytb_records, out_path, "fasta")

    print(f"Wrote {len(cytb_records)} sequences to {out_path}")
    for rec in cytb_records:
        print(f"{rec.id:12s} length = {len(rec.seq)} bp")

    return cytb_records


def main():
    records = write_cytb_fasta(ACCESSIONS, OUT_FASTA)
    print(records)
    breakpoint()


if __name__ == "__main__":
    main()
