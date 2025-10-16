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

GVCLASS_PREFIX = {
    "domain": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}


def _resolve_parquet_glob(parquet_glob: str | None) -> str:
    return parquet_glob or SETTINGS.parquet_glob


def _limit_clause(limit: int | None) -> str:
    return f"LIMIT {limit}" if limit is not None else ""


@lru_cache(maxsize=64)
def _cached_parquet(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


def _cache_frame(name: str) -> pd.DataFrame | None:
    if not SETTINGS.cache_enabled:
        return None
    path = SETTINGS.cache_dir / name
    if not path.exists():
        return None
    return _cached_parquet(str(path))


@lru_cache(maxsize=1)
def _metadata_frame() -> pd.DataFrame | None:
    path = SETTINGS.metadata_path
    if path is None:
        return None
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=["NA", "N/A"])
    df.columns = [col.lower().replace(" ", "_").replace("/", "_") for col in df.columns]
    if "genome" in df.columns:
        df = df.rename(columns={"genome": "dataset_id"})
    if "order_completeness" in df.columns:
        df["order_completeness"] = pd.to_numeric(df["order_completeness"], errors="coerce")
    return df


def _gvclass_rank_expr(level: str, fallback: str = "'Unassigned'") -> str:
    raw = _gvclass_rank_raw(level)
    return "COALESCE(" f"NULLIF(NULLIF({raw}, ''), 'NA'), {fallback})"


def _gvclass_rank_raw(level: str) -> str:
    prefix = GVCLASS_PREFIX[level]
    pattern = f"{prefix}_([^;]+)"
    return f"regexp_extract(taxonomy_majority, '{pattern}', 1)"


@contextmanager
def get_connection(
    parquet_glob: str | None = None,
) -> Iterator[catalog.duckdb.DuckDBPyConnection]:
    """Context manager that yields a read-only DuckDB connection."""

    resolved_glob = _resolve_parquet_glob(parquet_glob)
    conn = catalog.connect(db_path=":memory:", parquet_glob=resolved_glob, read_only=False)
    try:
        yield conn
    finally:
        catalog.close(conn)


def _frame_from_query(query: str, parquet_glob: str | None = None, **params) -> pd.DataFrame:
    with get_connection(parquet_glob) as conn:
        return conn.execute(query, params).df()


@lru_cache(maxsize=1)
def fetch_overview_metrics(
    parquet_glob: str | None = None,
    *,
    use_cache: bool = True,
) -> dict[str, float]:
    """Return high-level summary metrics used on the overview page."""

    if use_cache:
        cache = _cache_frame("overview.parquet")
        if cache is not None and not cache.empty:
            row = cache.iloc[0]
            return {
                "total_sequences": float(row["total_sequences"]),
                "unique_genomes": float(row["unique_genomes"]),
                "taxonomy_labels": float(row["taxonomy_labels"]),
                "gvclass_species": float(row["gvclass_species"]),
                "ani_clusters": float(row["ani_clusters"]),
                "avg_gc_nt": float(row["avg_gc_nt"]),
            }

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
    limit: int | None = 25,
    *,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Return taxonomy counts for the requested level and source."""

    cache_name = f"taxonomy_{source}_{level}.parquet"
    if use_cache:
        cache = _cache_frame(cache_name)
        if cache is not None:
            return cache.head(limit) if limit else cache.copy()

    if source == "gvclass":
        raw = _gvclass_rank_raw(level)
        limit_sql = _limit_clause(limit)
        df = _frame_from_query(
            f"""
            WITH base AS (
                SELECT
                    dataset_id,
                    {raw} AS raw_label
                FROM sequences
                WHERE seq_type = 'NT'
            )
        SELECT
            COALESCE(NULLIF(NULLIF(raw_label, ''), 'NA'), 'Unassigned') AS label,
            COUNT(DISTINCT dataset_id) AS genomes
        FROM base
        GROUP BY 1
        ORDER BY genomes DESC
        {limit_sql}
        """,
            parquet_glob,
        )
    else:
        column = PHYLO_COLUMNS.get(level)
        if column is None:
            raise ValueError(f"Unsupported taxonomy level: {level}")

        limit_sql = _limit_clause(limit)
        df = _frame_from_query(
            f"""
            WITH base AS (
                SELECT
                    dataset_id,
                    NULLIF({column}, '') AS raw_label
                FROM sequences
                WHERE seq_type = 'NT'
            )
        SELECT
            COALESCE(raw_label, 'Unassigned') AS label,
            COUNT(DISTINCT dataset_id) AS genomes
        FROM base
        GROUP BY 1
        ORDER BY genomes DESC
        {limit_sql}
        """,
            parquet_glob,
        )

    return df.head(limit) if limit else df


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
    limit: int | None = 30,
    *,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Return environment counts for the requested dimension."""

    column = ENVIRONMENT_COLUMNS.get(dimension)
    if column is None:
        raise ValueError(f"Unsupported environment dimension: {dimension}")

    cache_name = f"environment_{dimension}.parquet"
    if use_cache:
        cache = _cache_frame(cache_name)
        if cache is not None:
            return cache.head(limit) if limit else cache.copy()

    limit_sql = _limit_clause(limit)
    df = _frame_from_query(
        f"""
        WITH base AS (
            SELECT
                dataset_id,
                {column} AS raw_label
            FROM sequences
        )
        SELECT
            COALESCE(NULLIF(raw_label, ''), 'Unassigned') AS label,
            COUNT(DISTINCT dataset_id) AS genomes,
            COUNT(*) AS sequences
        FROM base
        GROUP BY 1
        ORDER BY genomes DESC
        {limit_sql}
        """,
        parquet_glob,
    )

    return df.head(limit) if limit else df


def fetch_genome_statistics(
    parquet_glob: str | None = None,
    *,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Return summary statistics for genomes with parsed taxonomy levels."""

    if use_cache:
        cache = _cache_frame("genome_statistics.parquet")
        if cache is not None:
            return cache.copy()

    df = _frame_from_query(
        """
        SELECT
            dataset_id,
            MAX(lenbp) AS genome_length,
            MAX(genecount) AS gene_count,
            AVG(gcperc) AS gc_percent,
            AVG(codingperc) AS coding_percent,
            MAX(contigs) AS contig_count,
            MAX(order_completeness) AS completeness,
            MAX(order_dup) AS contamination,
            MAX(taxonomy_majority) AS taxonomy_majority,
            MAX(ecosystem) AS ecosystem
        FROM sequences
        WHERE seq_type = 'NT'
        GROUP BY dataset_id
        """,
        parquet_glob,
    )

    if df.empty:
        return df

    metadata_df = _metadata_frame()
    if metadata_df is not None and "dataset_id" in metadata_df.columns:
        meta_subset = metadata_df[["dataset_id", "order_completeness"]].copy()
        df = df.merge(meta_subset, on="dataset_id", how="left", suffixes=("", "_meta"))
        if "order_completeness" in df.columns:
            df["completeness"] = df["order_completeness"].fillna(df["completeness"])
            df.drop(columns=["order_completeness"], inplace=True)
    df.rename(columns={"contamination": "order_dup"}, inplace=True)
    df["order_dup"] = df["order_dup"].astype(float)
    df.rename(columns={"order_dup": "contamination"}, inplace=True)

    taxonomy_map = {
        "taxonomy_domain": "d",
        "taxonomy_phylum": "p",
        "taxonomy_class": "c",
        "taxonomy_order": "o",
        "taxonomy_family": "f",
        "taxonomy_genus": "g",
        "taxonomy_species": "s",
    }

    for column, prefix in taxonomy_map.items():
        df[column] = (
            df["taxonomy_majority"]
            .str.extract(rf"{prefix}_([^;]+)")
            .iloc[:, 0]
            .fillna("Unassigned")
        )

    df["genome_length_mb"] = df["genome_length"] / 1e6
    df["gc_percent"] = df["gc_percent"].astype(float)
    df["coding_percent"] = df["coding_percent"].astype(float)
    df["completeness"] = df["completeness"].astype(float).fillna(0.0)
    df["contamination"] = df["contamination"].astype(float).fillna(0.0)
    df["contig_count"] = df["contig_count"].fillna(0).astype(int)

    return df


def fetch_annotation_matrix(
    field: Literal[
        "emapper_COG_category",
        "emapper_KEGG_Pathway",
        "emapper_PFAMs",
        "emapper_Description",
    ] = "emapper_COG_category",
    parquet_glob: str | None = None,
    limit: int | None = 25,
    *,
    level: Literal["order", "family", "genus", "class", "phylum", "domain"] = "order",
    use_cache: bool = True,
) -> pd.DataFrame:
    """Return annotation counts grouped by GVClass order for heatmap visualisations."""

    if level not in GVCLASS_PREFIX:
        raise ValueError(f"Unsupported GVClass level: {level}")

    cache_name = f"annotations_{field}_{level}.parquet"
    if use_cache:
        cache = _cache_frame(cache_name)
        if cache is not None:
            return cache.head(limit) if limit else cache.copy()

    raw_rank = _gvclass_rank_raw(level)

    limit_sql = _limit_clause(limit)
    df = _frame_from_query(
        f"""
        WITH base AS (
            SELECT
                {raw_rank} AS raw_order,
                NULLIF({field}, '') AS raw_annotation
            FROM sequences
            WHERE seq_type = 'AA'
        )
        SELECT
            COALESCE(NULLIF(raw_order, 'NA'), 'Unassigned') AS gvclass_order,
            COALESCE(raw_annotation, 'Unannotated') AS annotation_value,
            COUNT(*) AS sequences
        FROM base
        GROUP BY 1, 2
        ORDER BY sequences DESC
        {limit_sql}
        """,
        parquet_glob,
    )

    return df.head(limit) if limit else df


def fetch_cluster_summary(
    parquet_glob: str | None = None,
    limit: int | None = 50,
    *,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Summarise ANI clusters derived from skani results."""

    if use_cache:
        cache = _cache_frame("cluster_summary.parquet")
        if cache is not None:
            return cache.head(limit) if limit else cache.copy()

    raw_order = _gvclass_rank_raw("order")
    limit_sql = _limit_clause(limit)
    df = _frame_from_query(
        f"""
        WITH base AS (
            SELECT
                COALESCE(s_cluster, 'Unclustered') AS cluster_id,
                dataset_id,
                {raw_order} AS raw_order,
                is_skani_representative
            FROM sequences
            WHERE seq_type = 'NT'
        )
        SELECT
            cluster_id,
            COUNT(DISTINCT dataset_id) AS genomes,
            SUM(CASE WHEN is_skani_representative THEN 1 ELSE 0 END) AS representatives,
            COALESCE(MIN(NULLIF(raw_order, 'NA')), 'Unassigned') AS gvclass_order
        FROM base
        GROUP BY cluster_id
        ORDER BY genomes DESC
        {limit_sql}
        """,
        parquet_glob,
    )

    return df.head(limit) if limit else df
