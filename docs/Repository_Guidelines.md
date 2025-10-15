# Repository Guidelines

## Project Structure & Modules
All runtime code lives under `src/gvmagdb`, following a package-by-domain layout:
- `core/` manages DuckDB connectivity, catalog helpers, and shared utilities.
- `analyze/` hosts user-facing search and analysis commands (`run_skani.py`, `run_pyhmmer.py`, `tnf_distance.py`).
- `export/` includes the DIAMOND builder and FASTA exporters.
- `prepare/` wraps ingestion flows and annotation parsers.
Reusable fixtures for manual validation live in `tests/data`, while automated suites run from `tests/`. Large generated artifacts (`artifacts/` and `data/registry/`) are treated as read-only.

## Build, Test, and Dev Commands
```bash
pixi install                 # Resolve environments
pixi run install-package     # Editable install, exposes gvmagdb CLI
pixi run check               # fmt + lint + typecheck (Black, Ruff, mypy)
pixi run test                # Pytest suite
pixi run test-cov            # Coverage with HTML report
pixi run dashboard           # Launch Plotly Dash analytics UI (host/port configurable)
pixi run dashboard-cache     # Pre-compute dashboard aggregates (artifacts/analytics)
```
CLI tasks map 1:1 to modules: e.g., `pixi run build-diamond`, `pixi run export-faa`, `pixi run search-hmm`.

## Coding Style & Conventions
Formatting is enforced by Black (100-char lines) and Ruff; never hand-edit around those constraints. Use `snake_case` for modules and functions, `PascalCase` for classes, and short verb-based CLI command names (`diamond`, `extract`). Public functions should carry type hints; push business logic into `core` helpers to keep `cli.py` slim.

## Testing Expectations
Put unit tests under `tests/` with `test_*.py` names. Flag long runs with `@pytest.mark.slow` or `@pytest.mark.integration` to mirror existing markers. Prefer the lightweight data in `tests/data` when validating ingestion or search flows locally. Always run `pixi run test` before committing; use `pixi run test-cov` when code paths are hard to exercise.

## Commit & Pull Request Hygiene
Follow the `type: summary` pattern for commit subjects (`feat: add skani batch runner`). Each PR should:
1. Link issues or ticket IDs.
2. Describe behavior changes and affected CLI commands.
3. Include sample output or screenshots for user-visible changes.
4. State `pixi run check` and `pixi run test` results.

## Data & Security Notes
Avoid committing `artifacts/`, `data/registry/`, or any provenance-heavy downloads; document reproduction steps instead. Secrets and API credentials must come from environment variables (`.env` ignored by default) and never be stored in the repo. Update `docs/` when data schemas or CLI surfaces change.***
