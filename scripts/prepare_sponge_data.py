#!/usr/bin/env python3

"""
Download 18S rRNA sequences from key animal phyla and a choanoflagellate
outgroup, then align them to produce a FASTA file for phylogenetic analysis.

This dataset is inspired by Musser et al. (2021) "Profiling cellular diversity
in sponges informs animal cell type and nervous system evolution"
(Science, doi:10.1126/science.abj2949), which used the freshwater sponge
Spongilla lacustris to study the evolutionary origins of animal cell types and
the nervous system.  Building a phylogenetic tree with 18S rRNA sequences shows
where sponges sit in the animal tree of life, providing context for why they
are a key model organism for understanding early animal evolution.

Taxa included (one representative per major lineage):
    Monosiga_brevicollis   Choanoflagellata  – closest non-animal relatives (outgroup)
    Spongilla_lacustris    Porifera          – the study organism from Musser et al. 2021
    Trichoplax_sp          Placozoa          – simplest known animal (no tissues)
    Beroe_forskalii        Ctenophora        – comb jelly (early-branching)
    Nematostella_vectensis Cnidaria          – sea anemone (radially symmetric)
    Caenorhabditis_elegans Ecdysozoa         – nematode (bilaterian)
    Homo_sapiens           Vertebrata        – human (bilaterian, deuterostome)

The output is an aligned FASTA file compatible with scripts/run_clustering.py.

Dependencies:
    conda install -c conda-forge biopython mafft
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


# NCBI asks users to provide a contact email for Entrez queries.
Entrez.email = "constantin.pape@informatik.uni-goettingen.de"

# Verified NCBI 18S rRNA accession numbers for each taxon.
# All are complete or near-complete 18S sequences (~1,750–1,870 bp).
ACCESSIONS = {
    "Monosiga_brevicollis":   "AF100940.1",   # Choanoflagellata (outgroup)
    "Spongilla_lacustris":    "AF058945.1",   # Porifera – Musser et al. 2021
    "Trichoplax_sp":          "Z22783.1",     # Placozoa (1,787 bp)
    "Beroe_forskalii":        "AF293697.1",   # Ctenophora (1,802 bp)
    "Nematostella_vectensis": "AF254382.1",   # Cnidaria
    "Caenorhabditis_elegans": "NR_000054.1",  # Ecdysozoa – RefSeq (1,754 bp)
    "Homo_sapiens":           "NR_003286.4",  # Vertebrata – RefSeq (1,869 bp)
}

OUT_FASTA = Path("animal_18S_rRNA.fasta")
OUT_FASTA_ALIGNED = Path("animal_18S_rRNA_aligned.fasta")


def download_sequences(accessions: dict[str, str]) -> list[SeqRecord]:
    """Fetch sequences from NCBI GenBank by accession number."""
    ids = ",".join(accessions.values())
    print(f"Downloading {len(accessions)} sequences from NCBI …")

    with Entrez.efetch(db="nuccore", id=ids, rettype="fasta", retmode="text") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    # Re-label records with our short taxon names.
    # NCBI FASTA IDs are the bare accession, e.g. "AF100940.1"
    acc_to_label = {acc: label for label, acc in accessions.items()}
    relabelled = []
    for rec in records:
        rec_acc = rec.id  # e.g. "AF100940.1"
        label = acc_to_label.get(rec_acc)
        if label is None:
            # Fall back: match ignoring version suffix (e.g. "NR_000054" -> "NR_000054.1")
            base = rec_acc.rsplit(".", 1)[0]
            for key_acc, lbl in acc_to_label.items():
                if key_acc.rsplit(".", 1)[0] == base:
                    label = lbl
                    break
        if label is None:
            raise RuntimeError(f"Could not match downloaded record '{rec.id}' to a taxon label.")

        relabelled.append(
            SeqRecord(rec.seq, id=label, name=label, description="")
        )

    # Sort to match original order.
    order = list(accessions.keys())
    relabelled.sort(key=lambda r: order.index(r.id))
    return relabelled


def align_with_mafft(records: list[SeqRecord], aligned_path: Path) -> list[SeqRecord]:
    """Write an unaligned FASTA, run MAFFT, and return the aligned records."""
    if shutil.which("mafft") is None:
        raise RuntimeError(
            "MAFFT is not installed or not on PATH. "
            "Install it with:  conda install -c conda-forge mafft"
        )

    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as tmp:
        SeqIO.write(records, tmp, "fasta")
        tmp_path = Path(tmp.name)

    print("Aligning sequences with MAFFT …")
    result = subprocess.run(
        ["mafft", "--auto", "--quiet", str(tmp_path)],
        capture_output=True,
        text=True,
        check=True,
    )

    tmp_path.unlink()

    aligned_path.write_text(result.stdout)
    aligned_records = list(SeqIO.parse(aligned_path, "fasta"))
    print(f"Wrote aligned FASTA to {aligned_path}")
    return aligned_records


def main():
    records = download_sequences(ACCESSIONS)

    # Write unaligned sequences for reference.
    SeqIO.write(records, OUT_FASTA, "fasta")
    print(f"Wrote unaligned FASTA to {OUT_FASTA}")
    for rec in records:
        print(f"  {rec.id:30s}  {len(rec.seq):>5} bp")

    # Align and write the file expected by run_clustering.py.
    aligned = align_with_mafft(records, OUT_FASTA_ALIGNED)
    print(f"\n{len(aligned)} aligned sequences ready for phylogenetic analysis.")
    print(f"Pass '{OUT_FASTA_ALIGNED}' to run_clustering.py (update OUT_FASTA at the top).")


if __name__ == "__main__":
    main()
