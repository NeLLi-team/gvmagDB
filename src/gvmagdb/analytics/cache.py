"""Utilities to pre-compute analytics datasets for faster dashboards."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from gvmagdb.analytics import data_access
from gvmagdb.analytics.config import SETTINGS

LEVELS: tuple[str, ...] = ("domain", "phylum", "class", "order", "family", "genus", "species")
ENV_DIMENSIONS: tuple[str, ...] = (
    "ecosystem",
    "ecosystem_category",
    "ecosystem_type",
    "ecosystem_subtype",
    "habitat",
    "source",
)
ANNOTATION_FIELDS: tuple[str, ...] = (
    "emapper_COG_category",
    "emapper_KEGG_Pathway",
    "emapper_PFAMs",
    "emapper_Description",
)
ANNOTATION_LEVELS: tuple[str, ...] = ("order", "family", "genus", "class", "phylum", "domain")


def _write_parquet(name: str, df: pd.DataFrame) -> Path:
    SETTINGS.cache_dir.mkdir(parents=True, exist_ok=True)
    path = SETTINGS.cache_dir / f"{name}.parquet"
    df.to_parquet(path, index=False)
    return path


def build_cache(parquet_glob: str | None = None) -> None:
    """Compute and store cached analytics tables."""

    print("📦 Writing analytics cache to", SETTINGS.cache_dir)

    overview = pd.DataFrame([data_access.fetch_overview_metrics(parquet_glob, use_cache=False)])
    _write_parquet("overview", overview)

    taxonomy_frames: list[pd.DataFrame] = []
    for source in ("gvclass", "phylo"):
        for level in LEVELS:
            if source == "phylo" and level == "species":
                continue
            df = data_access.fetch_taxonomy_distribution(
                level,
                source,
                parquet_glob,
                limit=None,
                use_cache=False,
            )
            df.insert(0, "level", level)
            df.insert(0, "source", source)
            taxonomy_frames.append(df)
            _write_parquet(f"taxonomy_{source}_{level}", df)
    taxonomy_all = pd.concat(taxonomy_frames, ignore_index=True)
    _write_parquet("taxonomy", taxonomy_all)

    for dimension in ENV_DIMENSIONS:
        df = data_access.fetch_environment_distribution(
            dimension,
            parquet_glob,
            limit=None,
            use_cache=False,
        )
        _write_parquet(f"environment_{dimension}", df)

    genome_stats = data_access.fetch_genome_statistics(parquet_glob, use_cache=False)
    _write_parquet("genome_statistics", genome_stats)

    for field in ANNOTATION_FIELDS:
        level_frames: list[pd.DataFrame] = []
        for level in ANNOTATION_LEVELS:
            df = data_access.fetch_annotation_matrix(
                field,
                parquet_glob,
                limit=None,
                level=level,
                use_cache=False,
            )
            _write_parquet(f"annotations_{field}_{level}", df)
            df = df.copy()
            df.insert(0, "level", level)
            level_frames.append(df)
        combined = pd.concat(level_frames, ignore_index=True)
        _write_parquet(f"annotations_{field}", combined)

    cluster_summary = data_access.fetch_cluster_summary(parquet_glob, limit=1000)
    _write_parquet("cluster_summary", cluster_summary)

    print("✅ Analytics cache refreshed")


def main() -> None:
    build_cache()


if __name__ == "__main__":
    main()
