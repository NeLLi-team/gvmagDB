"""
TNF (tetranucleotide frequency) distance calculation.

Calculates TNF vectors for query sequences and compares them to database sequences
to find closest matches based on cosine distance or Euclidean distance.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from gvmagdb.core.catalog import connect


def calculate_tnf(sequence: str) -> np.ndarray:
    """
    Calculate tetranucleotide frequency (TNF) vector for a sequence.

    Args:
        sequence: DNA sequence string

    Returns:
        136-dimensional TNF vector (normalized frequencies)
    """
    sequence = sequence.upper()

    # Generate all 256 tetranucleotides (4^4)
    bases = ["A", "C", "G", "T"]
    tetranucleotides = []
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                for b4 in bases:
                    tetranucleotides.append(b1 + b2 + b3 + b4)

    # Count tetranucleotides
    tnf_counts = dict.fromkeys(tetranucleotides, 0)
    total = 0

    for i in range(len(sequence) - 3):
        tetra = sequence[i : i + 4]
        if all(base in bases for base in tetra):
            tnf_counts[tetra] += 1
            total += 1

    # Calculate frequencies
    tnf_freqs = np.array([tnf_counts[tn] / total if total > 0 else 0 for tn in tetranucleotides])

    # Use only 136 dimensions (removing palindromic redundancy)
    # Standard TNF-136 uses non-redundant tetranucleotides
    indices_136 = []
    seen = set()

    for i, tn in enumerate(tetranucleotides):
        rev_comp = "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[b] for b in reversed(tn))
        if tn not in seen and rev_comp not in seen:
            indices_136.append(i)
            seen.add(tn)
            seen.add(rev_comp)

    tnf_136 = tnf_freqs[indices_136[:136]]

    return tnf_136


def calculate_tnf_from_fasta(fasta_path: Path) -> tuple[str, np.ndarray]:
    """
    Calculate TNF vector from FASTA file (first sequence).

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Tuple of (sequence_id, tnf_vector)
    """
    with open(fasta_path) as f:
        record = next(SeqIO.parse(f, "fasta"))
        seq_id = record.id
        tnf = calculate_tnf(str(record.seq))

    return seq_id, tnf


def cosine_distance(v1: np.ndarray, v2: np.ndarray) -> float:
    """Calculate cosine distance between two vectors."""
    dot = np.dot(v1, v2)
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)

    if norm1 == 0 or norm2 == 0:
        return 1.0

    similarity = dot / (norm1 * norm2)
    distance = 1.0 - similarity

    return distance


def euclidean_distance(v1: np.ndarray, v2: np.ndarray) -> float:
    """Calculate Euclidean distance between two vectors."""
    return np.linalg.norm(v1 - v2)


def find_closest_by_tnf(
    query_fasta: Path,
    distance_metric: str = "cosine",
    top_k: int = 10,
    dataset_filter: str | None = None,
) -> pd.DataFrame:
    """
    Find closest genomes by TNF distance.

    Args:
        query_fasta: Path to query genome FASTA
        distance_metric: 'cosine' or 'euclidean'
        top_k: Number of top matches to return
        dataset_filter: Optional dataset_id to filter (for testing)

    Returns:
        DataFrame with closest matches and their distances
    """
    print("\n🧬 Calculating TNF distance...")
    print(f"   Query: {query_fasta}")
    print(f"   Distance metric: {distance_metric}")

    # Calculate query TNF
    query_id, query_tnf = calculate_tnf_from_fasta(query_fasta)
    print(f"   Query ID: {query_id}")
    print(f"   TNF dimensions: {len(query_tnf)}")

    # Connect to database
    conn = connect()

    try:
        # Query ONE representative contig per genome (longest with TNF)
        # Use genome-level length (lenbp) not contig length
        query = """
            WITH ranked_contigs AS (
                SELECT
                    dataset_id,
                    gid,
                    tnf_136,
                    length as contig_length,
                    lenbp as genome_length,
                    gcperc as genome_gc,
                    common_name as genome,
                    domain,
                    phylum,
                    class,
                    "order",
                    family,
                    genus,
                    gvclass_phylum,
                    gvclass_class,
                    gvclass_order,
                    gvclass_family,
                    gvclass_genus,
                    taxonomy_majority,
                    s_cluster,
                    skani_rep,
                    ROW_NUMBER() OVER (PARTITION BY dataset_id ORDER BY length DESC) as rn
                FROM sequences
                WHERE seq_type = 'NT' AND tnf_136 IS NOT NULL
        """

        if dataset_filter:
            query += f" AND dataset_id = '{dataset_filter}'"

        query += """
            )
            SELECT * FROM ranked_contigs WHERE rn = 1
        """

        print("\n📊 Querying database (one representative contig per genome)...")
        db_seqs = conn.execute(query).df()
        print(f"   Found {len(db_seqs)} genomes with TNF vectors")

        # Calculate distances
        distances = []

        for idx, row in db_seqs.iterrows():
            tnf_db = (
                json.loads(row["tnf_136"]) if isinstance(row["tnf_136"], str) else row["tnf_136"]
            )
            tnf_db = np.array(tnf_db)

            if distance_metric == "cosine":
                dist = cosine_distance(query_tnf, tnf_db)
            elif distance_metric == "euclidean":
                dist = euclidean_distance(query_tnf, tnf_db)
            else:
                raise ValueError(f"Unknown distance metric: {distance_metric}")

            distances.append(dist)

            if (idx + 1) % 1000 == 0:
                print(f"   Processed {idx + 1}/{len(db_seqs)} sequences...")

        db_seqs["tnf_distance"] = distances

        # Sort by distance and get top k
        results = db_seqs.sort_values("tnf_distance").head(top_k).copy()

        # Drop tnf_136 column for cleaner output
        results = results.drop(columns=["tnf_136"])

        print(f"\n✓ Found {len(results)} closest matches")

        return results

    finally:
        conn.close()


def compare_tnf_to_skani(
    query_fasta: Path,
    skani_results: pd.DataFrame,
    tnf_results: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compare TNF-based distances with skani ANI results.

    Args:
        query_fasta: Path to query genome
        skani_results: DataFrame from skani search
        tnf_results: DataFrame from TNF distance search

    Returns:
        Combined DataFrame showing both metrics
    """
    # Merge on dataset_id or gid
    if "Ref_name" in skani_results.columns:
        skani_results = skani_results.rename(columns={"Ref_name": "dataset_id"})

    combined = tnf_results.merge(
        skani_results[["dataset_id", "ANI", "Align_fraction_query"]], on="dataset_id", how="left"
    )

    # Sort by TNF distance
    combined = combined.sort_values("tnf_distance")

    return combined
