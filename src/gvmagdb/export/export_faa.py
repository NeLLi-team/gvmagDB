"""
Protein FASTA export module.

Exports protein sequences from Parquet/DuckDB with >{gid} headers
for use with DIAMOND, HMMER, and other search tools.
"""

from pathlib import Path

import duckdb

from gvmagdb.core.catalog import connect


def export_aa_fasta(
    output_faa: Path,
    dataset_id: str | None = None,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: Path | None = None,
    parquet_glob: str | None = None,
) -> int:
    """
    Export protein sequences to FASTA.

    Args:
        output_faa: Output FASTA file path
        dataset_id: Optional filter by dataset_id
        conn: Existing DuckDB connection (or will create one)
        db_path: DuckDB catalog path (if conn not provided)
        parquet_glob: Parquet glob pattern (if conn not provided)

    Returns:
        Number of sequences exported

    Notes:
        - FASTA headers are EXACTLY >{gid} with no spaces
        - This ensures trivial joins back to Parquet via gid
    """
    # Create connection if not provided
    close_conn = False
    if conn is None:
        conn = connect(
            db_path=db_path or Path("artifacts/duckdb/gvmagdb.duckdb"),
            parquet_glob=parquet_glob or "artifacts/parquet/sequences/**/*.parquet",
        )
        close_conn = True

    try:
        # Query protein sequences
        query = "SELECT gid, sequence FROM sequences WHERE seq_type = 'AA'"
        params = []

        if dataset_id:
            query += " AND dataset_id = ?"
            params.append(dataset_id)

        query += " ORDER BY gid"

        # Ensure output directory exists
        output_faa.parent.mkdir(parents=True, exist_ok=True)

        # Write FASTA
        count = 0
        with open(output_faa, "w") as f:
            for row in conn.execute(query, params).fetchall():
                gid, sequence = row
                f.write(f">{gid}\n")
                # Write sequence in 60-char lines
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i : i + 60] + "\n")
                count += 1

        print(f"✓ Exported {count} protein sequences to {output_faa}")
        return count

    finally:
        if close_conn:
            conn.close()


def export_aa_by_genomes(
    output_dir: Path,
    genome_ids: list[str],
    conn: duckdb.DuckDBPyConnection | None = None,
) -> dict[str, int]:
    """
    Export protein sequences per genome.

    Args:
        output_dir: Directory for per-genome FASTA files
        genome_ids: List of genome IDs to export
        conn: Optional DuckDB connection

    Returns:
        Dict mapping genome_id -> count of sequences
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    counts = {}

    for genome_id in genome_ids:
        output_faa = output_dir / f"{genome_id}.faa"
        count = export_aa_fasta(output_faa, dataset_id=genome_id, conn=conn)
        counts[genome_id] = count

    return counts
