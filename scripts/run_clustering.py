from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt


OUT_DIST = Path("primate_mt_cytb_pdist.csv")
OUT_TREE = Path("primate_mt_cytb_upgma.png")
OUT_FASTA = Path("primate_mt_cytb.fasta")


def p_distance(seq1, seq2):
    """
    Simple p-distance between two aligned sequences:
    fraction of compared positions that differ.

    Positions with ambiguous bases or gaps are ignored.
    """
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()

    if len(seq1) != len(seq2):
        raise ValueError(
            f"Sequences must have equal length, got {len(seq1)} and {len(seq2)}. "
            "For this CYTB dataset they should normally already match."
        )

    valid = {"A", "C", "G", "T"}

    n_compared = 0
    n_different = 0

    for a, b in zip(seq1, seq2):
        if a in valid and b in valid:
            n_compared += 1
            n_different += int(a != b)

    if n_compared == 0:
        raise ValueError("No comparable positions found.")

    return n_different / n_compared


def compute_distance_matrix(records):
    labels = [rec.id for rec in records]
    n = len(records)

    distances = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            d = p_distance(records[i].seq, records[j].seq)
            distances[i, j] = d
            distances[j, i] = d

    return pd.DataFrame(distances, index=labels, columns=labels)


def plot_upgma(distance_df, out_path):
    """
    Plot an UPGMA-like tree using average linkage.

    Note:
        SciPy's method='average' is average-linkage hierarchical clustering,
        which corresponds to the same clustering idea used in UPGMA.
    """
    condensed = squareform(distance_df.values)
    Z = linkage(condensed, method="average")

    plt.figure(figsize=(7, 4))
    dendrogram(Z, labels=list(distance_df.index))
    plt.ylabel("p-distance")
    plt.title("Primate mitochondrial CYTB: UPGMA / average-linkage tree")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    print(f"Wrote tree plot to {out_path}")


def main():
    records = list(SeqIO.parse(OUT_FASTA, "fasta"))
    distance_df = compute_distance_matrix(records)

    print("\nPairwise p-distance matrix:")
    print(distance_df.round(4))

    plot_upgma(distance_df, OUT_TREE)


if __name__ == "__main__":
    main()
