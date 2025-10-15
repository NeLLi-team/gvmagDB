"""
Data access helpers for the analytics dashboards.

These utilities provide small, memoised wrappers around DuckDB queries so that
Plotly Dash callbacks can remain lean and focussed on presentation logic.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from functools import lru_cache
from typing import Literal

import pandas as pd

from gvmagdb.core import catalog

from .config import SETTINGS

GVCLASS_COLUMNS = {
    "domain": "gvclass_domain",
    "phylum": "gvclass_phylum",
    "class": "gvclass_class",
    "order": "gvclass_order",
    "family": "gvclass_family",
    "genus": "gvclass_genus",
    "species": "gvclass_species",
}

PHYLO_COLUMNS = {
    "domain": "domain",
    "phylum": "phylum",
    "class": "class",
    "order": '"order"',
    "family": "family",
    "genus": "genus",
}

ENVIRONMENT_COLUMNS = {
    "ecosystem": "ecosystem",
    "ecosystem_category": "ecosystem_category",
    "ecosystem_type": "ecosystem_type",
    "ecosystem_subtype": "ecosystem_subtype",
    "habitat": "habitat",
    "source": "source",
}


def _resolve_parquet_glob(parquet_glob: str | None) -> str:
    return parquet_glob or SETTINGS.parquet_glob


@contextmanager
def get_connection(
    parquet_glob: str | None = None,
) -> Iterator[catalog.duckdb.DuckDBPyConnection]:
    """Context manager that yields a read-only DuckDB connection."""

    conn = catalog.connect(parquet_glob=_resolve_parquet_glob(parquet_glob), read_only=False)
    try:
        yield conn
    finally:
        catalog.close(conn)


def _frame_from_query(query: str, parquet_glob: str | None = None, **params) -> pd.DataFrame:
    with get_connection(parquet_glob) as conn:
        return conn.execute(query, params).df()


@lru_cache(maxsize=1)
def fetch_overview_metrics(parquet_glob: str | None = None) -> dict[str, float]:
    """Return high-level summary metrics used on the overview page."""

    row = _frame_from_query(
        """
        SELECT
            COUNT(*) AS total_sequences,
            COUNT(DISTINCT dataset_id) AS unique_genomes,
            COUNT(DISTINCT taxonomy_majority) AS taxonomy_labels,
            COUNT(DISTINCT gvclass_species) AS gvclass_species,
            COUNT(DISTINCT s_cluster) AS ani_clusters,
            AVG(gc_content) FILTER (WHERE seq_type = 'NT' AND gc_content IS NOT NULL) AS avg_gc_nt
        FROM sequences
        """,
        parquet_glob,
    ).iloc[0]

    return {
        "total_sequences": float(row["total_sequences"]),
        "unique_genomes": float(row["unique_genomes"]),
        "taxonomy_labels": float(row["taxonomy_labels"]),
        "gvclass_species": float(row["gvclass_species"]),
        "ani_clusters": float(row["ani_clusters"]),
        "avg_gc_nt": float(row["avg_gc_nt"]) if row["avg_gc_nt"] is not None else 0.0,
    }


def fetch_taxonomy_distribution(
    level: Literal["domain", "phylum", "class", "order", "family", "genus", "species"],
    source: Literal["gvclass", "phylo"] = "gvclass",
    parquet_glob: str | None = None,
    limit: int = 25,
) -> pd.DataFrame:
    """Return taxonomy counts for the requested level and source."""

    column_map = GVCLASS_COLUMNS if source == "gvclass" else PHYLO_COLUMNS
    column = column_map.get(level)
    if column is None:
        raise ValueError(f"Unsupported taxonomy level: {level}")

    return _frame_from_query(
        f"""
        SELECT
            COALESCE({column}, 'Unknown') AS label,
            COUNT(DISTINCT dataset_id) AS genomes
        FROM sequences
        WHERE seq_type = 'NT'
        GROUP BY label
        ORDER BY genomes DESC
        LIMIT {limit}
        """,
        parquet_glob,
    )


def fetch_environment_distribution(
    dimension: Literal[
        "ecosystem",
        "ecosystem_category",
        "ecosystem_type",
        "ecosystem_subtype",
        "habitat",
        "source",
    ],
    parquet_glob: str | None = None,
    limit: int = 30,
) -> pd.DataFrame:
    """Return environment counts for the requested dimension."""

    column = ENVIRONMENT_COLUMNS.get(dimension)
    if column is None:
        raise ValueError(f"Unsupported environment dimension: {dimension}")

    return _frame_from_query(
        f"""
        SELECT
            COALESCE({column}, 'Unknown') AS label,
            COUNT(DISTINCT dataset_id) AS genomes,
            COUNT(*) AS sequences
        FROM sequences
        GROUP BY label
        ORDER BY genomes DESC
        LIMIT {limit}
        """,
        parquet_glob,
    )


def fetch_genome_statistics(parquet_glob: str | None = None) -> pd.DataFrame:
    """Return summary statistics for genome length, gene count, GC%, coding density."""

    return _frame_from_query(
        """
        SELECT
            dataset_id,
            MAX(lenbp) AS genome_length,
            MAX(genecount) AS gene_count,
            AVG(gcperc) AS gc_percent,
            AVG(codingperc) AS coding_percent,
            MAX(taxonomy_majority) AS taxonomy_majority,
            MAX(ecosystem) AS ecosystem,
            MAX(order_completeness) AS order_completeness
        FROM sequences
        WHERE seq_type = 'NT'
        GROUP BY dataset_id
        """,
        parquet_glob,
    )


def fetch_annotation_matrix(
    field: Literal[
        "emapper_COG_category",
        "emapper_KEGG_Pathway",
        "emapper_PFAMs",
        "emapper_Description",
    ] = "emapper_COG_category",
    parquet_glob: str | None = None,
    limit: int = 25,
) -> pd.DataFrame:
    """Return annotation counts grouped by GVClass order for heatmap visualisations."""

    return _frame_from_query(
        f"""
        SELECT
            COALESCE(gvclass_order, 'Unknown') AS gvclass_order,
            COALESCE({field}, 'Unannotated') AS annotation_value,
            COUNT(*) AS sequences
        FROM sequences
        WHERE seq_type = 'AA'
        GROUP BY gvclass_order, annotation_value
        ORDER BY sequences DESC
        LIMIT {limit}
        """,
        parquet_glob,
    )


def fetch_cluster_summary(parquet_glob: str | None = None, limit: int = 50) -> pd.DataFrame:
    """Summarise ANI clusters derived from skani results."""

    return _frame_from_query(
        f"""
        SELECT
            COALESCE(s_cluster, 'Unclustered') AS cluster_id,
            COUNT(DISTINCT dataset_id) AS genomes,
            SUM(CASE WHEN is_skani_representative THEN 1 ELSE 0 END) AS representatives,
            MAX(gvclass_order) AS gvclass_order
        FROM sequences
        WHERE seq_type = 'NT'
        GROUP BY cluster_id
        ORDER BY genomes DESC
        LIMIT {limit}
        """,
        parquet_glob,
    )
