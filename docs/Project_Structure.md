# Project Structure

```
gvmagDB/
├── docs/                  # Contributor & user documentation (this folder)
├── src/gvmagdb/           # Installable Python package
│   ├── cli.py             # Entry point wired to `gvmagdb` console script
│   ├── analyze/           # Search & analytics commands
│   ├── core/              # Shared infrastructure (DuckDB catalog, TNF math, registries)
│   ├── export/            # DIAMOND builder and FASTA/TSV exporters
│   └── prepare/           # Ingestion flows and annotation parsing
├── tests/                 # Pytest suites (unit + integration)
│   └── data/              # Lightweight fixtures used by tests and demos
├── artifacts/             # Generated Parquet + search indices (ignored)
│   ├── parquet/           # Primary columnar storage
│   ├── products/          # DIAMOND, FASTA, skani references
│   └── duckdb/            # Local DuckDB catalogs (optional)
├── data/registry/         # Product registry cache (generated)
├── bin/                   # `gvmagdb` launcher for system PATH usage
├── pixi.toml              # Environment + task runner definitions
├── pyproject.toml         # Packaging, dependencies, tooling configuration
└── legacy/                # Archived scripts and previous implementations
```

## Module Overview
- `analyze/analyze_genome.py`: combines TNF and ANI comparisons for genome placement.
- `analyze/run_*`: wrappers for DIAMOND, PyHMMER, and Skani search flows.
- `core/catalog.py`: central DuckDB connection factory and schema helpers.
- `core/tnf.py`: Numba-accelerated TNF feature generation.
- `core/registry.py`: tracks generated product versions and checksums.
- `export/`: CLI-accessible utilities for building DIAMOND databases and exporting FASTA outputs.
- `prepare/ingest.py`: sequential ingestion pipeline (FASTA → Parquet).
- `prepare/ingest_parallel.py`: MPI-style high-throughput ingestion harness.
- `prepare/parse_annotations.py`: ProteinOrtho / eggNOG ingestion helpers.

## Data & Legacy Notes
Everything in `artifacts/` and `data/registry/` is derived output and intentionally excluded from version control; regenerate them via the CLI workflows documented in `docs/Database_Workflows.md`. Anything under `legacy/` is retained for historical reference and should not receive new changes unless you are porting functionality into the supported code paths. Update `docs/` whenever new modules or directories are introduced so the high-level overview stays current.***
