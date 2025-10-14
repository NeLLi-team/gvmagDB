# Database Schema Reference

## Storage Layout
All canonical data lands in Parquet files under `artifacts/parquet/`, partitioned by `dataset_id`. DuckDB is used as the query engine; `core/catalog.py` loads the Parquet directory as an external table set, so no manual imports are required.

```
artifacts/parquet/
├── sequences/          # Sequence-level metadata (AA + NT)
├── genomes/            # Genome-level annotations and quality metrics
├── taxonomy/           # GVClass and original taxonomy assignments
├── orthogroups/        # ProteinOrtho clusters and related metrics
└── annotations/        # Functional annotations (eggNOG, KEGG, PFAM, GO)
```

## Key Tables
- **`genomes`** (72 columns): genome taxonomy snapshots, clustering labels (ANI, phylogenetic, MCL), quality metrics, environment metadata, IMG/JGI provenance, and summary statistics (length, GC%, gene counts).
- **`sequences`** (23 columns): per-sequence identifiers, feature lengths, TNF signatures, coding strand, translation tables, and representative genome flags.
- **`annotations`**: functional labels (`nog_category`, `kegg_module`, `pfam_domain`, `go_terms`, `ec_number`) mapped by `sequence_id`. Each record carries confidence scores and source metadata.
- **`orthogroups`**: ProteinOrtho-derived clusters with group-level statistics (size, representative, consensus function), keyed by `orthogroup_id`.
- **`taxonomy`**: original submissions (domain → genus) and GVClass re-assignments (strict and majority calls across seven ranks).

## Access Patterns
Typical DuckDB session:
```python
from gvmagdb.core.catalog import connect
conn = connect(read_only=True)
df = conn.execute(\"\"\"SELECT genome_name, ani_cluster, completeness
                    FROM genomes
                    WHERE completeness > 90
                    ORDER BY completeness DESC
                    LIMIT 50\"\"\").df()
```

Use `JOIN` against `sequences` for per-protein metadata, and `JOIN` against `annotations` for functional signal. When adding columns, update `core/catalog.py` so new fields are exposed in the default views and ensure ingestion scripts populate values for both sequential and parallel paths.***
