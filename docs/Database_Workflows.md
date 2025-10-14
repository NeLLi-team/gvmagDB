# Database Workflows

This guide walks through the four core operational flows for gvmagDB: provisioning the database, running local analytics, querying with new sequences from external environments, and planning hosted dashboards.

## 1. Provision the Database
1. **Install tooling**
   ```bash
   pixi install
   pixi run install-package    # exposes the `gvmagdb` CLI
   ```
2. **Download genomic assets** from the shared object store (newlineages bucket or internal mirror):
   ```bash
   mkdir -p ingestion_data
   curl -LO https://newlineages.com/data/GVMAGS_V2_data.tar.gz
   tar -xzf GVMAGS_V2_data.tar.gz -C ingestion_data
   ```
   The tarball unpacks `gvmagsV2all.fna`, `gvmagsV2all.faa`, `Updated_naming_Sept2025.tsv`, and auxiliary annotation folders.
3. **Ingest into Parquet**
   ```bash
   gvmagdb ingest-cmd \
     --fna ingestion_data/gvmagsV2all.fna \
     --faa ingestion_data/gvmagsV2all.faa \
     --metadata ingestion_data/Updated_naming_Sept2025.tsv \
     --output artifacts/parquet/sequences
   ```
4. **Verify catalog state**
   ```bash
   gvmagdb stats
   ```
   Expect ~18k genomes and 5M sequences. If the counts differ, re-run ingestion and confirm the tarball was fully extracted.

## 2. Local Analytics & Bulk Exports
- **Export every protein FASTA for the order *Imitervirales***:
  ```bash
  pixi run python - <<'PY'
from pathlib import Path
from gvmagdb.core.catalog import connect
from gvmagdb.export.export_faa import export_aa_by_genomes

conn = connect()
rows = conn.execute("""
    SELECT DISTINCT dataset_id
    FROM sequences
    WHERE "order" = 'Imitervirales'
""").fetchall()
ids = [row[0] for row in rows]
export_aa_by_genomes(Path("exports/imitervirales"), ids, conn=conn)
conn.close()
PY
  ```
  The command emits one `.faa` file per genome inside `exports/imitervirales/`.
- **Genome size distribution plot** (Plotly):
  ```bash
  pixi run python - <<'PY'
import plotly.express as px
from gvmagdb.core.catalog import connect

conn = connect()
df = conn.execute("""
    SELECT dataset_id, MAX(lenbp) AS genome_bp, "order"
    FROM sequences
    WHERE seq_type = 'NT'
    GROUP BY dataset_id, "order"
""").df()
fig = px.histogram(df, x="genome_bp", color="order",
                   nbins=60, title="Genome size distribution (bp)")
fig.write_html("exports/genome_size_distribution.html")
conn.close()
PY
  ```
  Open the HTML file locally to explore interactive histograms.

## 3. Remote Sequence Queries
1. **Pre-build search indices** on the workstation:
   ```bash
   gvmagdb export-fna-cmd --output artifacts/products/fna/gvmagdb.fna
   gvmagdb export-faa-cmd --output artifacts/products/faa/gvmagdb.faa
   gvmagdb build-diamond --output artifacts/products/dmnd/gvmagdb.dmnd --threads 32
   ```
2. **Sync the artifacts** (`rsync` or `aws s3 cp`) to the HPC or cloud instance where analyses will run.
3. **ANI clustering with skani** for new genomes:
   ```bash
   gvmagdb skani \
     --query my_genomes.fna \
     --ref /data/gvmagdb.fna \
     --ani-threshold 85 \
     --species-threshold 95 \
     --threads 32
   ```
   Inspect `*_vs_gvmagdb/gvmagdb_ani.tsv` and `cluster_membership.txt` to determine companion clusters.
4. **Protein homology with DIAMOND**:
   ```bash
   gvmagdb diamond \
     --query proteins.faa \
     --db /data/gvmagdb.dmnd \
     --evalue 1e-10 \
     --threads 32 \
     --extract
   ```
   Hit sequences land in `gvmagdb_hits.faa`; metadata-rich tables appear in `gvmagdb_hits.tsv`.
5. **Profile HMM search**:
   ```bash
   gvmagdb hmmsearch \
     --hmm markers.hmm \
     --target /data/gvmagdb.faa \
     --threads 32 \
     --extract
   ```
   Use `gvmagdb extract --hits gvmagdb_hits.tsv --fasta filtered_hits.faa --seq-type AA` to repackage results for downstream workflows.

## 4. Plotly Dashboards (Planned)
- **Goal**: deliver interactive analytics via Plotly Dash hosted on `newlineages.com` through Cloudflare Tunnel.
- **Proposed stack**:
  - Dash app in `src/gvmagdb/analytics/dashboard.py`, reading DuckDB views in read-only mode.
  - Pre-computed parquet summaries (`docs/releases/<version>.json`) for lightweight rendering.
  - Deployment using `pixi run python -m gvmagdb.analytics.dashboard --port 8050`.
  - Cloudflare Tunnel configuration stored under `infrastructure/cloudflare/`, mapping `https://newlineages.com/gvmagdb`.
- **Next steps**: define API boundaries (read-only endpoints), produce mock layouts (genome size histograms, taxonomy sunburst, search hit explorer), and add automated screenshots to release notes.

Keep this document updated as commands evolve or the dashboard implementation lands.***
