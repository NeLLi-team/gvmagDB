# GVMAG Genome Database

High-performance viral and microbial genome repository with Parquet storage, DuckDB querying, and integrated search tooling.

## Overview
- **Source layout**: installable package under `src/gvmagdb`, CLI entry point `gvmagdb`.
- **Generated data**: written to `artifacts/` (`parquet/`, `products/`, `duckdb/`) to keep the repo clean.
- **Documentation**: full contributor and operations guides live in `docs/` (start with `docs/README.md`).

## Quick Start
```bash
pixi install                     # create/solve environment
pixi run install-package         # editable install with console script
pixi run ingest \                # ingest genomes + metadata from FASTA/TSV
  --fna ingestion_data/gvmagsV2all.fna \
  --faa ingestion_data/gvmagsV2all.faa \
  --metadata ingestion_data/Updated_naming_Sept2025.tsv
pixi run dashboard-cache         # pre-compute dashboard summaries (optional but recommended)
pixi run stats                   # report high-level database metrics
pixi run dashboard               # launch the interactive analytics UI
```
Daily development helpers:
- `pixi run fmt`, `pixi run lint`, `pixi run typecheck` (or `pixi run check` for the full bundle).
- `pixi run build-diamond`, `pixi run export-fna`, `pixi run export-faa` to refresh search products.
- `pixi run search-diamond`, `search-hmm`, `skani` for query workflows.

## Database & Workflow Map
```mermaid
flowchart LR
    subgraph Source Data
        A[Tarball: GVMAGS_V2_data.tar.gz]
    end
    subgraph Ingestion
        B[CLI: gvmagdb ingest-cmd]
        C[Annotations Loader]
    end
    subgraph Storage
        D[artifacts/parquet<br/>DuckDB views]
        E[artifacts/products<br/>(FASTA, DMND, FNA)]
        F[artifacts/duckdb<br/>local catalog]
    end
    subgraph Workflows
        G[Query: gvmagdb diamond/hmmsearch/skani]
        H[Analytics: docs/Database_Workflows<br/>Plotly Dash app]
    end

    A --> B
    B --> C --> D
    B --> E
    D --> F
    D --> G
    E --> G
    G --> H
    D --> H
```

- Ingestion pulls FASTA/metadata from `ingestion_data/`, enriches with ProteinOrtho & eggNOG annotations, and writes Parquet partitions plus search products into `artifacts/`.
- Query commands read DuckDB views (no duplication) and reuse generated FASTA/DMND artifacts.
- The Plotly Dash dashboards consume the same artifacts via read-only DuckDB connections (see `docs/Database_Workflows.md` and `docs/Repository_Guidelines.md`).

## Analytics Cache

The dashboards scan ~5M sequences; precomputing aggregates keeps page loads responsive.

- Run `pixi run dashboard-cache` after ingestion or whenever the database changes. Results live under `artifacts/analytics/` (override via `GVMAGDB_ANALYTICS_CACHE_DIR`).
- At runtime the data-access layer automatically reads cached Parquet files. Set `GVMAGDB_ANALYTICS_CACHE=false` to bypass caching during development.

## Further Reading
- [Repository Guidelines](docs/Repository_Guidelines.md) – coding style, testing, PR hygiene.
- [Database Workflows](docs/Database_Workflows.md) – provisioning, analytics, remote querying, dashboard roadmap.
- [Database Schema Reference](docs/Database_Schema.md) – Parquet partitions and table descriptions.
- [Distribution Workflow](docs/Distribution_Workflow.md) – release checklist and artifact publishing.
