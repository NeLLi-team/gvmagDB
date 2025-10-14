"""
FASTA → Parquet ingestion module.

Writes sequences directly to partitioned Parquet files (PRIMARY storage).
Preserves ALL metadata columns from TSV.
"""

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO

from gvmagdb.core.tnf import calculate_tnf, tnf_to_list
from gvmagdb.prepare.parse_annotations import load_all_annotations

# Parquet schema with ALL metadata fields preserved
SEQUENCES_SCHEMA = pa.schema(
    [
        # Core sequence fields
        ("gid", pa.string()),  # dataset_id|header
        ("dataset_id", pa.string()),  # from genome column
        ("header", pa.string()),  # original FASTA ID
        ("seq_type", pa.string()),  # 'AA' or 'NT'
        ("length", pa.int32()),
        ("gc_content", pa.float32()),  # NT only
        ("sequence", pa.string()),
        ("tnf_136", pa.list_(pa.float32())),  # NT only, 136-dim
        # ALL metadata columns from Updated_naming_Sept2025.tsv
        ("genome", pa.string()),
        ("common_name", pa.string()),
        ("phylum", pa.string()),
        ("class", pa.string()),
        ("order", pa.string()),
        ("family", pa.string()),
        ("genus", pa.string()),
        ("s_cluster", pa.string()),
        ("skani_rep", pa.string()),
        ("pdmtree", pa.string()),
        ("pdm_cluster", pa.string()),
        ("mcl_rep", pa.string()),
        ("mcl_on_pdm46tree", pa.string()),
        ("order_completeness", pa.float64()),
        ("gvog4_unique", pa.int64()),
        ("gvog8_unique", pa.int64()),
        ("gvog8_total", pa.int64()),
        ("gvog8_dup", pa.float64()),
        ("mcp_total", pa.int64()),
        ("mirus_unique", pa.int64()),
        ("mirus_total", pa.int64()),
        ("mirus_dup", pa.float64()),
        ("mrya_unique", pa.int64()),
        ("mrya_total", pa.int64()),
        ("lenbp", pa.int64()),
        ("gcperc", pa.float64()),
        ("genecount", pa.int64()),
        ("codingperc", pa.float64()),
        ("ttable", pa.string()),
        ("domain", pa.string()),
        ("sequencing_status", pa.string()),
        ("study_name", pa.string()),
        ("genome_name___sample_name", pa.string()),
        ("sequencing_center", pa.string()),
        ("img_genome_id", pa.string()),
        ("gold_analysis_project_type", pa.string()),
        ("assembly_method", pa.string()),
        ("pubmed_id", pa.string()),
        ("ecosystem", pa.string()),
        ("ecosystem_category", pa.string()),
        ("ecosystem_subtype", pa.string()),
        ("ecosystem_type", pa.string()),
        ("specific_ecosystem", pa.string()),
        ("habitat", pa.string()),
        ("source", pa.string()),
        # GVClass taxonomy assignments (prefixed with gvclass_)
        ("taxonomy_majority", pa.string()),
        ("taxonomy_strict", pa.string()),
        ("gvclass_species", pa.string()),
        ("gvclass_genus", pa.string()),
        ("gvclass_family", pa.string()),
        ("gvclass_order", pa.string()),
        ("gvclass_class", pa.string()),
        ("gvclass_phylum", pa.string()),
        ("gvclass_domain", pa.string()),
        # GVClass additional stats
        ("avgdist", pa.float64()),
        ("order_dup", pa.float64()),
        ("order_weighted_completeness", pa.float64()),
        ("order_confidence_score", pa.float64()),
        ("phage_unique", pa.int64()),
        ("phage_total", pa.int64()),
        ("cellular_unique", pa.int64()),
        ("cellular_total", pa.int64()),
        ("cellular_dup", pa.float64()),
        ("contigs", pa.int64()),
        # Protein annotations (only for skani representative genomes)
        ("is_skani_representative", pa.bool_()),
        ("orthogroup_id", pa.string()),
        ("is_singleton", pa.bool_()),
        # eggNOG-mapper functional annotations
        ("emapper_seed_ortholog", pa.string()),
        ("emapper_evalue", pa.float64()),
        ("emapper_score", pa.float64()),
        ("emapper_eggNOG_OGs", pa.string()),
        ("emapper_max_annot_lvl", pa.string()),
        ("emapper_COG_category", pa.string()),
        ("emapper_Description", pa.string()),
        ("emapper_Preferred_name", pa.string()),
        ("emapper_GOs", pa.string()),
        ("emapper_EC", pa.string()),
        ("emapper_KEGG_ko", pa.string()),
        ("emapper_KEGG_Pathway", pa.string()),
        ("emapper_KEGG_Module", pa.string()),
        ("emapper_KEGG_Reaction", pa.string()),
        ("emapper_KEGG_rclass", pa.string()),
        ("emapper_BRITE", pa.string()),
        ("emapper_KEGG_TC", pa.string()),
        ("emapper_CAZy", pa.string()),
        ("emapper_BiGG_Reaction", pa.string()),
        ("emapper_PFAMs", pa.string()),
    ]
)


def parse_metadata_tsv(metadata_path: Path) -> tuple[dict, dict]:
    """
    Parse metadata TSV and merge with gvclass results.

    Args:
        metadata_path: Path to Updated_naming_Sept2025.tsv

    Returns:
        Tuple of (metadata_dict, skani_rep_dict)
        - metadata_dict: genome_id -> metadata dict (including gvclass data)
        - skani_rep_dict: genome_id -> is_representative (boolean)
    """
    # Read base metadata
    df = pd.read_csv(metadata_path, sep="\t", keep_default_na=False, na_values=["NA", "N/A"])

    # Normalize column names (lowercase, replace spaces with underscores)
    df.columns = [col.lower().replace(" ", "_").replace("/", "_") for col in df.columns]

    # Convert numeric columns explicitly with proper types
    int_cols = [
        "gvog4_unique",
        "gvog8_unique",
        "gvog8_total",
        "mcp_total",
        "mirus_unique",
        "mirus_total",
        "mrya_unique",
        "mrya_total",
        "lenbp",
        "genecount",
    ]
    float_cols = ["order_completeness", "gvog8_dup", "mirus_dup", "gcperc", "codingperc"]

    for col in int_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").round().astype("Int64")

    for col in float_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Convert all non-numeric columns to strings to avoid type mismatches
    all_numeric_cols = int_cols + float_cols
    string_cols = [col for col in df.columns if col not in all_numeric_cols]
    for col in string_cols:
        df[col] = df[col].astype(str)

    # Load gvclass data
    gvclass_path = Path("ingestion_data/gvclass/all_results.tsv")
    if gvclass_path.exists():
        print(f"📊 Loading gvclass data from {gvclass_path}...")
        gvclass_df = pd.read_csv(
            gvclass_path, sep="\t", keep_default_na=False, na_values=["NA", "N/A"]
        )

        # Rename columns: query -> genome, taxonomy columns with gvclass_ prefix
        gvclass_df = gvclass_df.rename(columns={"query": "genome"})

        # Prefix taxonomy columns with gvclass_
        taxonomy_rename = {
            "species": "gvclass_species",
            "genus": "gvclass_genus",
            "family": "gvclass_family",
            "order": "gvclass_order",
            "class": "gvclass_class",
            "phylum": "gvclass_phylum",
            "domain": "gvclass_domain",
        }
        gvclass_df = gvclass_df.rename(columns=taxonomy_rename)

        # Convert gvclass numeric columns with proper types
        gvclass_int_cols = [
            "gvog4_unique",
            "gvog8_unique",
            "gvog8_total",
            "mcp_total",
            "mirus_unique",
            "mirus_total",
            "mrya_unique",
            "mrya_total",
            "phage_unique",
            "phage_total",
            "cellular_unique",
            "cellular_total",
            "contigs",
            "lenbp",
            "genecount",
        ]
        gvclass_float_cols = [
            "avgdist",
            "order_dup",
            "order_completeness",
            "order_weighted_completeness",
            "order_confidence_score",
            "gvog8_dup",
            "mirus_dup",
            "cellular_dup",
            "gcperc",
            "codingperc",
        ]

        for col in gvclass_int_cols:
            if col in gvclass_df.columns:
                gvclass_df[col] = (
                    pd.to_numeric(gvclass_df[col], errors="coerce").round().astype("Int64")
                )

        for col in gvclass_float_cols:
            if col in gvclass_df.columns:
                gvclass_df[col] = pd.to_numeric(gvclass_df[col], errors="coerce")

        # Convert gvclass string columns
        gvclass_all_numeric = gvclass_int_cols + gvclass_float_cols
        gvclass_string_cols = [
            col for col in gvclass_df.columns if col not in gvclass_all_numeric and col != "genome"
        ]
        for col in gvclass_string_cols:
            gvclass_df[col] = gvclass_df[col].astype(str)

        # Merge on genome ID, using gvclass values for overlapping columns
        # First, drop overlapping columns from base df (except genome)
        overlap_cols = [
            "order_completeness",
            "gvog4_unique",
            "gvog8_unique",
            "gvog8_total",
            "gvog8_dup",
            "mcp_total",
            "mirus_unique",
            "mirus_total",
            "mirus_dup",
            "mrya_unique",
            "mrya_total",
            "lenbp",
            "gcperc",
            "genecount",
            "codingperc",
            "ttable",
        ]
        cols_to_drop = [
            col for col in overlap_cols if col in df.columns and col in gvclass_df.columns
        ]
        df = df.drop(columns=cols_to_drop)

        # Merge, keeping all genomes from base metadata
        df = df.merge(gvclass_df, on="genome", how="left", suffixes=("", "_gvclass"))

        print(f"✓ Merged gvclass data for {len(gvclass_df)} genomes")
    else:
        print(f"⚠ Warning: gvclass data not found at {gvclass_path}, skipping")
        # Add empty gvclass columns
        gvclass_cols = [
            "taxonomy_majority",
            "taxonomy_strict",
            "gvclass_species",
            "gvclass_genus",
            "gvclass_family",
            "gvclass_order",
            "gvclass_class",
            "gvclass_phylum",
            "gvclass_domain",
            "avgdist",
            "order_dup",
            "order_weighted_completeness",
            "order_confidence_score",
            "phage_unique",
            "phage_total",
            "cellular_unique",
            "cellular_total",
            "cellular_dup",
            "contigs",
        ]
        for col in gvclass_cols:
            df[col] = None

    # Build lookup dicts
    metadata_dict = {}
    skani_rep_dict = {}

    for _, row in df.iterrows():
        genome_id = row["genome"]
        metadata_dict[genome_id] = row.to_dict()

        # Determine if genome is skani representative
        # skani_rep column contains the representative genome ID
        # If genome == skani_rep, then this genome IS the representative
        skani_rep_id = str(row.get("skani_rep", ""))
        is_rep = genome_id == skani_rep_id
        skani_rep_dict[genome_id] = is_rep

    return metadata_dict, skani_rep_dict


def detect_seq_type(sequence: str) -> str:
    """Detect if sequence is NT or AA."""
    seq_upper = sequence.upper()
    nt_bases = set("ACGTUN")
    unique_chars = set(seq_upper)

    # If >95% of unique chars are ACGTUN, it's nucleotide
    nt_count = sum(1 for c in unique_chars if c in nt_bases)
    if nt_count / max(1, len(unique_chars)) > 0.95:
        return "NT"
    return "AA"


def parse_fasta_header(header: str) -> tuple[str, str]:
    """
    Parse FASTA header to extract dataset_id and local header.

    Supports formats:
    1. genome_id|contig_id
    2. contig_id (extract genome_id from _NUMBER pattern)

    Returns:
        (dataset_id, header)
    """
    if "|" in header:
        dataset_id, local_header = header.split("|", 1)
        return dataset_id, local_header
    else:
        # Try to extract genome_id from pattern: genome_id_NUMBER
        parts = header.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0], header
        else:
            return header, header


def calculate_gc(sequence: str) -> float:
    """Calculate GC content percentage."""
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    total = len(seq_upper)
    if total == 0:
        return 0.0
    return (gc_count / total) * 100.0


def fasta_to_parquet(
    fasta_path: Path,
    metadata_dict: dict[str, dict[str, object]],
    skani_rep_dict: dict[str, bool],
    parquet_root: Path,
    orthogroup_dict: dict[str, str] | None = None,
    singleton_set: set[str] | None = None,
    eggnog_dict: dict[str, dict[str, object | None]] | None = None,
    batch_size: int = 50000,
) -> None:
    """
    Ingest FASTA file to partitioned Parquet.

    Args:
        fasta_path: Path to FNA or FAA file
        metadata_dict: Genome metadata dictionary
        skani_rep_dict: Skani representative status dictionary
        parquet_root: Root directory for Parquet files
        orthogroup_dict: Protein → orthogroup_id mapping (optional)
        singleton_set: Set of singleton protein IDs (optional)
        eggnog_dict: Protein → eggNOG annotations mapping (optional)
        batch_size: Records per batch for memory efficiency
    """
    # Define eggNOG annotation columns
    emapper_columns = [
        "seed_ortholog",
        "evalue",
        "score",
        "eggNOG_OGs",
        "max_annot_lvl",
        "COG_category",
        "Description",
        "Preferred_name",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    ]

    print(f"\n📖 Parsing FASTA: {fasta_path}...")
    rows: list[dict[str, object]] = []
    skipped = 0
    processed = 0

    def flush_batch():
        if not rows:
            return
        # Group by dataset_id to avoid partition limit (max 1024 partitions per write)
        df = pd.DataFrame(rows)
        for _dataset_id, group_df in df.groupby("dataset_id"):
            table = pa.Table.from_pandas(group_df, schema=SEQUENCES_SCHEMA)
            pq.write_to_dataset(
                table,
                root_path=str(parquet_root),
                partition_cols=["dataset_id"],
                compression="zstd",
                compression_level=6,
                use_dictionary=["dataset_id", "seq_type"],
            )
        rows.clear()

    for record in SeqIO.parse(fasta_path, "fasta"):
        dataset_id, local_header = parse_fasta_header(record.id)

        # Skip if no metadata for this genome
        if dataset_id not in metadata_dict:
            skipped += 1
            continue

        sequence = str(record.seq).upper().replace("U", "T")
        seq_type = detect_seq_type(sequence)
        gid = f"{dataset_id}|{local_header}"

        # Calculate NT-specific fields
        gc_content = None
        tnf = None
        if seq_type == "NT":
            gc_content = calculate_gc(sequence)
            tnf_array = calculate_tnf(sequence)
            tnf = tnf_to_list(tnf_array)

        # Determine if genome is skani representative
        is_rep = skani_rep_dict.get(dataset_id, False)

        # Initialize annotation fields
        orthogroup_id = None
        is_singleton = None
        emapper_annotations: dict[str, object | None] = {}

        # Lookup annotations if genome is representative
        if is_rep:
            # Check orthogroup membership
            if orthogroup_dict and gid in orthogroup_dict:
                orthogroup_id = orthogroup_dict[gid]
                is_singleton = False
            elif singleton_set and gid in singleton_set:
                orthogroup_id = None
                is_singleton = True
            # else: protein not in proteinortho results (is_singleton stays None)

            # Check eggNOG annotations
            if eggnog_dict and gid in eggnog_dict:
                emapper_annotations = eggnog_dict[gid]

        # Build row with ALL metadata and annotations
        meta = metadata_dict[dataset_id]
        row = {
            "gid": gid,
            "dataset_id": dataset_id,
            "header": local_header,
            "seq_type": seq_type,
            "length": len(sequence),
            "gc_content": gc_content,
            "sequence": sequence,
            "tnf_136": tnf,
            # All metadata fields (types already correct from parse_metadata_tsv)
            **meta,
            # Annotation fields
            "is_skani_representative": is_rep,
            "orthogroup_id": orthogroup_id,
            "is_singleton": is_singleton,
        }

        # Add eggNOG annotation fields
        for col in emapper_columns:
            row[f"emapper_{col}"] = emapper_annotations.get(col, None)

        rows.append(row)
        processed += 1

        if len(rows) >= batch_size:
            print(f"  Writing batch ({processed} sequences)...")
            flush_batch()

    # Flush remaining
    if rows:
        print(f"  Writing final batch ({processed} sequences)...")
        flush_batch()

    print("\n✓ Ingestion complete:")
    print(f"  - Processed: {processed} sequences")
    print(f"  - Skipped (no metadata): {skipped}")
    print(f"  - Output: {parquet_root}")


def ingest_all(
    metadata_path: Path,
    fna_path: Path,
    faa_path: Path,
    parquet_root: Path = Path("artifacts/parquet/sequences"),
) -> None:
    """
    Ingest both nucleotide and protein sequences with annotations.

    Args:
        metadata_path: Path to metadata TSV
        fna_path: Path to nucleotide FASTA
        faa_path: Path to protein FASTA
        parquet_root: Root for Parquet output
    """
    parquet_root.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("GVMAG Database Ingestion")
    print("=" * 70)

    # Load metadata
    print(f"\n📖 Loading metadata from {metadata_path}...")
    metadata_dict, skani_rep_dict = parse_metadata_tsv(metadata_path)
    print(f"✓ Loaded metadata for {len(metadata_dict)} genomes")
    print(f"  - Skani representatives: {sum(skani_rep_dict.values())}")

    # Load annotations (for representative genomes only)
    print("\n📊 Loading protein annotations...")
    proteinortho_tsv = Path("ingestion_data/proteinortho/gvmags_db_trial4.proteinortho.tsv")
    singletons_tsv = Path("ingestion_data/proteinortho/gvmags_db_trial4.singletons.tsv")
    eggnog_tsv = Path("ingestion_data/annotation/emapper_out.emapper.annotations")

    orthogroup_dict: dict[str, str] | None = None
    singleton_set: set[str] | None = None
    eggnog_dict: dict[str, dict[str, object | None]] | None = None

    if all([proteinortho_tsv.exists(), singletons_tsv.exists(), eggnog_tsv.exists()]):
        orthogroup_dict, singleton_set, eggnog_dict = load_all_annotations(
            proteinortho_tsv, singletons_tsv, eggnog_tsv
        )
        print("✓ Annotations loaded successfully")
    else:
        print("⚠ Warning: Annotation files not found, skipping annotation integration")
        if not proteinortho_tsv.exists():
            print(f"  - Missing: {proteinortho_tsv}")
        if not singletons_tsv.exists():
            print(f"  - Missing: {singletons_tsv}")
        if not eggnog_tsv.exists():
            print(f"  - Missing: {eggnog_tsv}")

    # Ingest nucleotides
    print("\n🧬 Ingesting nucleotide sequences...")
    fasta_to_parquet(
        fna_path,
        metadata_dict,
        skani_rep_dict,
        parquet_root,
        orthogroup_dict,
        singleton_set,
        eggnog_dict,
    )

    # Ingest proteins
    print("\n🧬 Ingesting protein sequences...")
    fasta_to_parquet(
        faa_path,
        metadata_dict,
        skani_rep_dict,
        parquet_root,
        orthogroup_dict,
        singleton_set,
        eggnog_dict,
    )

    print("\n" + "=" * 70)
    print("✓ Ingestion complete!")
    print("=" * 70)
