"""
Extract sequences and metadata from search results.

Takes search results (TSV with gid column) and extracts sequences
with enriched metadata from Parquet.
"""

import json
from pathlib import Path

import pandas as pd

from gvmagdb.core.catalog import connect


def extract_sequences(
    hits_tsv: Path,
    output_fasta: Path | None = None,
    output_tsv: Path | None = None,
    output_json: Path | None = None,
    gid_column: str = "sseqid",
    seq_type: str | None = None,
) -> pd.DataFrame:
    """
    Extract sequences and metadata from search results.

    Args:
        hits_tsv: TSV file with search results containing gid column
        output_fasta: Output FASTA file (optional)
        output_tsv: Output TSV file with full metadata (optional)
        output_json: Output JSON file with full metadata (optional)
        gid_column: Column name containing gids (default: 'sseqid' for BLAST/DIAMOND)
        seq_type: Filter by sequence type ('AA' or 'NT', None = both)

    Returns:
        DataFrame with extracted sequences and metadata
    """
    print(f"\n📥 Extracting sequences from {hits_tsv}...")

    # Read hits
    hits = pd.read_csv(hits_tsv, sep="\t")

    if gid_column not in hits.columns:
        # Try alternative column names
        if "target_name" in hits.columns:
            gid_column = "target_name"
        elif "gid" in hits.columns:
            gid_column = "gid"
        else:
            raise ValueError(
                f"Could not find gid column. Available columns: {hits.columns.tolist()}"
            )

    gids = hits[gid_column].unique().tolist()
    print(f"   Found {len(gids)} unique gids")

    # Connect to catalog
    conn = connect()

    try:
        # Build query
        placeholders = ",".join(["?"] * len(gids))
        seq_filter = f"AND seq_type = '{seq_type}'" if seq_type else ""

        query = f"""
            SELECT
                gid,
                dataset_id,
                header,
                seq_type,
                length,
                gc_content,
                sequence,
                genome,
                common_name,
                phylum,
                class,
                order,
                family,
                genus,
                s_cluster,
                ecosystem,
                ecosystem_category,
                ecosystem_type,
                habitat,
                lenbp,
                gcperc,
                genecount
            FROM sequences
            WHERE gid IN ({placeholders}) {seq_filter}
        """

        results = conn.execute(query, gids).df()
        print(f"✓ Extracted {len(results)} sequences")

        # Output FASTA
        if output_fasta:
            output_fasta.parent.mkdir(parents=True, exist_ok=True)
            with open(output_fasta, "w") as f:
                for _, row in results.iterrows():
                    # Enriched FASTA header with metadata
                    header_parts = [
                        f">{row['gid']}",
                        f"genome={row['genome']}",
                        f"phylum={row['phylum']}",
                        f"length={row['length']}",
                    ]
                    f.write(" | ".join(header_parts) + "\n")

                    # Write sequence in 60-char lines
                    seq = row["sequence"]
                    for i in range(0, len(seq), 60):
                        f.write(seq[i : i + 60] + "\n")

            print(f"✓ FASTA written to {output_fasta}")

        # Output TSV
        if output_tsv:
            output_tsv.parent.mkdir(parents=True, exist_ok=True)
            results.to_csv(output_tsv, sep="\t", index=False)
            print(f"✓ TSV written to {output_tsv}")

        # Output JSON
        if output_json:
            output_json.parent.mkdir(parents=True, exist_ok=True)
            with open(output_json, "w") as f:
                json.dump(results.to_dict(orient="records"), f, indent=2)
            print(f"✓ JSON written to {output_json}")

        return results

    finally:
        conn.close()


def extract_by_gids(
    gids: list[str],
    output_fasta: Path | None = None,
    seq_type: str | None = None,
) -> pd.DataFrame:
    """
    Extract sequences directly by list of gids.

    Args:
        gids: List of gids to extract
        output_fasta: Output FASTA file (optional)
        seq_type: Filter by sequence type ('AA' or 'NT', None = both)

    Returns:
        DataFrame with extracted sequences and metadata
    """
    print(f"\n📥 Extracting {len(gids)} sequences...")

    # Connect to catalog
    conn = connect()

    try:
        # Build query
        placeholders = ",".join(["?"] * len(gids))
        seq_filter = f"AND seq_type = '{seq_type}'" if seq_type else ""

        query = f"""
            SELECT
                gid,
                dataset_id,
                seq_type,
                length,
                sequence,
                genome,
                phylum
            FROM sequences
            WHERE gid IN ({placeholders}) {seq_filter}
        """

        results = conn.execute(query, gids).df()
        print(f"✓ Extracted {len(results)} sequences")

        # Output FASTA if requested
        if output_fasta:
            output_fasta.parent.mkdir(parents=True, exist_ok=True)
            with open(output_fasta, "w") as f:
                for _, row in results.iterrows():
                    f.write(f">{row['gid']} | genome={row['genome']} | phylum={row['phylum']}\n")
                    seq = row["sequence"]
                    for i in range(0, len(seq), 60):
                        f.write(seq[i : i + 60] + "\n")

            print(f"✓ FASTA written to {output_fasta}")

        return results

    finally:
        conn.close()
