"""
DIAMOND BLASTP search with join-back to Parquet.

Searches query proteins against DIAMOND database and joins results
back to sequence metadata via gid.
"""

import subprocess
from pathlib import Path

import pandas as pd

from gvmagdb.core.catalog import connect

# Standard BLAST output format 6 columns
BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def run_blastp(
    query_faa: Path,
    dmnd_db: Path,
    output_tsv: Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 100,
    threads: int = 8,
) -> pd.DataFrame:
    """
    Run DIAMOND BLASTP search.

    Args:
        query_faa: Query protein FASTA
        dmnd_db: DIAMOND database path (.dmnd file)
        output_tsv: Output TSV path
        evalue: E-value threshold
        max_target_seqs: Max hits per query
        threads: Number of threads

    Returns:
        DataFrame with BLAST results
    """
    print("\n🔍 Running DIAMOND BLASTP...")
    print(f"   Query: {query_faa}")
    print(f"   Database: {dmnd_db}")
    print(f"   E-value: {evalue}")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Remove .dmnd extension if present (diamond expects base name)
    dmnd_base = str(dmnd_db.with_suffix("")) if dmnd_db.suffix == ".dmnd" else str(dmnd_db)

    cmd = [
        "diamond",
        "blastp",
        "--query",
        str(query_faa),
        "--db",
        dmnd_base,
        "--out",
        str(output_tsv),
        "--outfmt",
        "6",
        "--evalue",
        str(evalue),
        "--max-target-seqs",
        str(max_target_seqs),
        "--threads",
        str(threads),
    ]

    subprocess.run(cmd, check=True)

    # Parse results
    if output_tsv.exists() and output_tsv.stat().st_size > 0:
        df = pd.read_csv(output_tsv, sep="\t", names=BLAST_COLUMNS)
        print(f"✓ Found {len(df)} hits")
        return df
    else:
        print("✓ No hits found")
        return pd.DataFrame(columns=BLAST_COLUMNS)


def blastp_with_metadata(
    query_faa: Path,
    dmnd_db: Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 100,
    threads: int = 8,
) -> pd.DataFrame:
    """
    Run DIAMOND BLASTP and join with sequence metadata.

    Args:
        query_faa: Query protein FASTA
        dmnd_db: DIAMOND database path
        evalue: E-value threshold
        max_target_seqs: Max hits per query
        threads: Number of threads

    Returns:
        DataFrame with BLAST results + metadata
    """
    # Run BLAST
    output_tsv = Path("results/diamond_results.tsv")
    results = run_blastp(query_faa, dmnd_db, output_tsv, evalue, max_target_seqs, threads)

    if results.empty:
        return results

    print("\n🔗 Joining results with metadata...")

    # Connect to catalog
    conn = connect()

    try:
        # Get unique target gids (sseqid IS gid because we exported with >{gid})
        target_gids = results["sseqid"].unique().tolist()

        # Build placeholders for SQL IN clause
        placeholders = ",".join(["?"] * len(target_gids))

        # Query metadata for all target sequences
        query = f"""
            SELECT
                gid,
                dataset_id,
                length,
                genome,
                phylum,
                class,
                order,
                family,
                genus,
                ecosystem,
                habitat
            FROM sequences
            WHERE gid IN ({placeholders})
        """

        metadata = conn.execute(query, target_gids).df()

        # Join results with metadata
        annotated = results.merge(metadata, left_on="sseqid", right_on="gid", how="left")

        print(f"✓ Annotated {len(annotated)} hits with metadata")
        return annotated

    finally:
        conn.close()
