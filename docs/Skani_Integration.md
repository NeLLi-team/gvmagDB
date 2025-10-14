# Skani Integration Notes

Skani provides fast approximate ANI calculations and clustering support for gvmagDB. The CLI entry point is exposed via `gvmagdb analyze skani` (implemented in `src/gvmagdb/analyze/run_skani.py`).

## Runtime Workflow
1. Ensure the reference database is registered (`artifacts/parquet/sequences/` must contain representative genomes with TNF signatures).
2. Run:
   ```bash
   gvmagdb skani --query my_genomes.fna --threads 16 --min-ani 0.85
   ```
3. The command:
   - Streams query genomes into Skani.
   - Joins results with DuckDB metadata for taxonomy and environment fields.
   - Emits TSV and JSON summaries under `<query>_vs_gvmagdb/`.

## Configuration Flags
- `--min-ani`: discard hits below the threshold (default 0.80).
- `--max-results`: limit matches per query (default 100).
- `--metadata-cols`: comma-separated columns to attach (defaults to `genome_name`, `gvclass_taxon`, `completeness`).
- `--output`: override the output directory.

## Extensibility Guidelines
- All Skani-specific logic lives in `run_skani.py`; keep CLI parsing thin and encapsulate data joins in helper functions.
- Shared constants (e.g., metadata column defaults) belong in `core/registry.py` to avoid duplication with other search tools.
- When Skani upgrades introduce new flags, map them through `click.option` definitions and document usage here.

## Testing
- Lightweight regression tests should mock Skani subprocess calls; integration tests can be gated behind `@pytest.mark.integration` and depend on the small FASTA fixtures in `tests/data`.
- Use `pixi run test -m "not integration"` for quick checks, then run full integration suites before release.

Update this document whenever output schemas change or new metadata joins are introduced, so downstream users can adjust parsing scripts accordingly.***
