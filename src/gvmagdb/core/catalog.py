"""
DuckDB catalog with views over Parquet files.

This module creates a thin DuckDB catalog that queries Parquet directly
without ingesting data. Supports local and remote (S3/GCS) Parquet files.
"""

from pathlib import Path

import duckdb
import pandas as pd


def connect(
    db_path: Path | str | None = Path("artifacts/duckdb/gvmagdb.duckdb"),
    parquet_glob: str = "artifacts/parquet/sequences/**/*.parquet",
    s3_config: dict | None = None,
    read_only: bool = False,
) -> duckdb.DuckDBPyConnection:
    """
    Create DuckDB connection with view over Parquet files.

    Args:
        db_path: Path to DuckDB catalog file
        parquet_glob: Glob pattern for Parquet files
        s3_config: Optional S3 configuration (e.g., {'s3_access_key_id': '...', ...})
        read_only: Open database in read-only mode

    Returns:
        DuckDB connection with 'sequences' view configured

    Example:
        >>> conn = connect()
        >>> result = conn.execute("SELECT COUNT(*) FROM sequences WHERE seq_type='NT'").fetchone()
        >>> print(f"Total NT sequences: {result[0]}")
    """
    # Ensure db directory exists
    if db_path is None or str(db_path) == ":memory:":
        conn = duckdb.connect(database=":memory:")
    else:
        db_path = Path(db_path)
        db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = duckdb.connect(str(db_path), read_only=read_only)

    # Configure for performance
    conn.execute("SET memory_limit='16GB'")
    conn.execute("SET threads=8")
    conn.execute("SET max_memory='16GB'")

    # Install and load httpfs for S3/GCS support
    conn.execute("INSTALL httpfs")
    conn.execute("LOAD httpfs")

    # Configure S3 if provided
    if s3_config:
        for key, value in s3_config.items():
            conn.execute(f"SET {key}='{value}'")

    # Create view over Parquet files
    # This does NOT ingest data - it's a view that queries Parquet directly
    conn.execute(
        f"""
        CREATE OR REPLACE VIEW sequences AS
        SELECT * FROM parquet_scan('{parquet_glob}')
    """
    )

    return conn


def close(conn: duckdb.DuckDBPyConnection) -> None:
    """Close DuckDB connection."""
    if conn:
        conn.close()


def get_stats(conn: duckdb.DuckDBPyConnection) -> dict:
    """
    Get database statistics.

    Args:
        conn: DuckDB connection

    Returns:
        Dictionary with statistics
    """
    stats = {}

    # Total sequences
    result = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
    stats["total_sequences"] = result[0]

    # By seq_type
    result = conn.execute(
        """
        SELECT seq_type, COUNT(*) as count
        FROM sequences
        GROUP BY seq_type
    """
    ).fetchall()
    stats["by_seq_type"] = {row[0]: row[1] for row in result}

    # By dataset
    result = conn.execute(
        """
        SELECT dataset_id, COUNT(*) as count
        FROM sequences
        GROUP BY dataset_id
        ORDER BY count DESC
        LIMIT 10
    """
    ).fetchall()
    stats["top_10_datasets"] = {row[0]: row[1] for row in result}

    # Unique genomes
    result = conn.execute("SELECT COUNT(DISTINCT dataset_id) FROM sequences").fetchone()
    stats["unique_genomes"] = result[0]

    # Total sequence length
    result = conn.execute("SELECT SUM(length) FROM sequences").fetchone()
    stats["total_bases"] = result[0]

    # Average GC content (NT only)
    result = conn.execute(
        """
        SELECT AVG(gc_content)
        FROM sequences
        WHERE seq_type='NT' AND gc_content IS NOT NULL
    """
    ).fetchone()
    stats["avg_gc_content"] = result[0] if result[0] else 0.0

    return stats


def query_by_gid(conn: duckdb.DuckDBPyConnection, gid: str) -> dict | None:
    """
    Retrieve sequence by global ID (gid).

    Args:
        conn: DuckDB connection
        gid: Global sequence ID (dataset_id|header)

    Returns:
        Dictionary with sequence data, or None if not found
    """
    result = conn.execute("SELECT * FROM sequences WHERE gid = ?", [gid]).fetchone()

    if not result:
        return None

    # Get column names
    columns = [desc[0] for desc in conn.description]
    return dict(zip(columns, result, strict=False))


def query_by_dataset(
    conn: duckdb.DuckDBPyConnection,
    dataset_id: str,
    seq_type: str | None = None,
    limit: int | None = None,
) -> list[dict]:
    """
    Query sequences by dataset_id.

    Args:
        conn: DuckDB connection
        dataset_id: Dataset ID to query
        seq_type: Optional filter by 'AA' or 'NT'
        limit: Optional limit on results

    Returns:
        List of sequence dictionaries
    """
    query = "SELECT * FROM sequences WHERE dataset_id = ?"
    params = [dataset_id]

    if seq_type:
        query += " AND seq_type = ?"
        params.append(seq_type)

    if limit:
        query += f" LIMIT {limit}"

    result = conn.execute(query, params).fetchall()
    columns = [desc[0] for desc in conn.description]

    return [dict(zip(columns, row, strict=False)) for row in result]


def export_sequences_to_df(
    conn: duckdb.DuckDBPyConnection,
    seq_type: str,
    dataset_id: str | None = None,
) -> pd.DataFrame:
    """
    Export sequences to pandas DataFrame.

    Args:
        conn: DuckDB connection
        seq_type: 'AA' or 'NT'
        dataset_id: Optional filter by dataset

    Returns:
        Pandas DataFrame
    """
    query = "SELECT * FROM sequences WHERE seq_type = ?"
    params = [seq_type]

    if dataset_id:
        query += " AND dataset_id = ?"
        params.append(dataset_id)

    return conn.execute(query, params).df()
