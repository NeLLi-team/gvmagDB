# Distribution Workflow

## Packaging Targets
- **Editable installs** (`pip install -e .` via `pixi run install-package`) for local development.
- **CLI distribution** through the `gvmagdb` console script defined in `pyproject.toml`.
- **Data bundles** (Parquet + DIAMOND/FASTA products) produced outside of Git history and shipped as separate artifacts (tarballs, object storage).

## Release Checklist
1. Run `pixi run check` and `pixi run test-cov`; ensure coverage does not regress.
2. Regenerate search indices as needed:
   ```bash
   pixi run build-diamond
   pixi run export-faa
   pixi run export-fna
   ```
3. Validate schema integrity:
   ```bash
   pixi run python -m gvmagdb.prepare.ingest --dry-run
   pixi run python -m gvmagdb.core.catalog --check
   ```
4. Update changelog content in `docs/` and confirm command help text (`gvmagdb --help`) reflects the release.
5. Bump `version` in `pyproject.toml` and `pixi.toml` where necessary.

## Artifact Publishing
- **Python Package**: create an sdist/wheel with `python -m build`, then upload via `twine upload dist/*`.
- **Data Products**: bundle the `artifacts/parquet/` partitions and `artifacts/products/` assets using `tar` or `pixi run python -m gvmagdb.export.bundle`. Upload to the designated object store (S3/GCS). Record checksums in `docs/releases/<version>.json`.
- **Container Images (optional)**: build from the repository root using `docker build -t gvmagdb:<tag> .` and push to the shared registry.

## Post-Release Tasks
- Tag the commit: `git tag -a vX.Y.Z -m "gvmagDB vX.Y.Z"` followed by `git push origin vX.Y.Z`.
- Announce the release notes with a summary of CLI changes, schema updates, and data refresh dates.
- Archive old artifacts and link their locations in `docs/releases/README.md` for reproducibility.

Maintain this workflow document as steps evolve, especially when onboarding new storage locations or changing the packaging toolchain.***
