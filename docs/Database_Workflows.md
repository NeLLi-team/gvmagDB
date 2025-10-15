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

## 4. Plotly Dashboards
- **Pre-compute analytics cache** *(recommended)*:
  ```bash
  pixi run dashboard-cache
  ```
  Cached Parquet files are written to `artifacts/analytics/` (override via `GVMAGDB_ANALYTICS_CACHE_DIR`).
- **Launch locally**:
  ```bash
  pixi run dashboard
  ```
  The app runs at `http://127.0.0.1:8050` by default and automatically loads assets from `dashboard/assets/`.
- **Pages included**:
  1. Overview – KPI cards plus quick taxonomy/environment snapshots.
  2. Taxonomy Explorer – compare GVClass vs phylogenomic assignments, drill down by rank.
  3. Environment & Ecosystems – visualise genome counts by ecosystem, habitat, and source.
  4. Genome Statistics – interactively inspect genome length, GC%, coding density, and completeness.
  5. Annotations – heatmaps for COG/KEGG/PFAM signals across GVClass orders.
  6. Clusters & Quality – summarise ANI clusters and representative coverage.
- **Data access**: the Dash app calls the shared DuckDB catalog (read-only) via `gvmagdb.analytics.data_access`, so no additional ingestion is required as long as `artifacts/parquet/` and `artifacts/products/` are present. Cached Parquet extracts are consumed automatically when available.
- **Deployment**: wrap the app with Gunicorn or Waitress (`python -m gvmagdb.analytics.app`) and expose via Tailscale Funnel (`tailscale funnel --bg --https=443 --set-path=/ 127.0.0.1:8050/`) or your preferred reverse proxy. Configure `GVMAGDB_PARQUET_GLOB`, `GVMAGDB_DASH_HOST`, `GVMAGDB_DASH_PORT`, and `GVMAGDB_DASH_DEBUG` environment variables as needed.
- **Enhancements**: enable authentication or Cloudflare Access before exposing the app publicly.

Keep this document updated as commands evolve or the dashboard implementation lands.***
