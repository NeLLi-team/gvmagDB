#!/usr/bin/env python
"""
Parallel FASTA → Parquet ingestion using multiprocessing.

Splits FASTA files into chunks and processes them in parallel using up to 16 CPUs.
"""
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO

from gvmagdb.core.tnf import calculate_tnf, tnf_to_list
from gvmagdb.prepare.ingest import (
    SEQUENCES_SCHEMA,
    calculate_gc,
    detect_seq_type,
    parse_fasta_header,
    parse_metadata_tsv,
)
from gvmagdb.prepare.parse_annotations import load_all_annotations


def process_fasta_chunk(args):
    """
    Process a chunk of FASTA records with annotations.

    Args:
        args: Tuple of (chunk_records, metadata_dict, skani_rep_dict, parquet_root,
                        orthogroup_dict, singleton_set, eggnog_dict, chunk_id)

    Returns:
        (processed_count, skipped_count)
    """
    (
        chunk_records,
        metadata_dict,
        skani_rep_dict,
        parquet_root,
        orthogroup_dict,
        singleton_set,
        eggnog_dict,
        chunk_id,
    ) = args

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

    rows: list[dict[str, object]] = []
    skipped = 0
    processed = 0
    batch_size = 10000

    print(f"  [Chunk {chunk_id}] Processing {len(chunk_records)} sequences...")

    for record in chunk_records:
        dataset_id, local_header = parse_fasta_header(record.id)

        # Skip if no metadata
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

            # Check eggNOG annotations
            if eggnog_dict and gid in eggnog_dict:
                emapper_annotations = eggnog_dict[gid]

        # Build row with all metadata and annotations
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
            **meta,
            "is_skani_representative": is_rep,
            "orthogroup_id": orthogroup_id,
            "is_singleton": is_singleton,
        }

        # Add eggNOG annotation fields
        for col in emapper_columns:
            row[f"emapper_{col}"] = emapper_annotations.get(col, None)

        rows.append(row)
        processed += 1

        # Write in batches
        if len(rows) >= batch_size:
            flush_batch(rows, parquet_root)
            rows.clear()

    # Flush remaining
    if rows:
        flush_batch(rows, parquet_root)

    print(f"  [Chunk {chunk_id}] ✓ Processed {processed}, skipped {skipped}")
    return processed, skipped


def flush_batch(rows, parquet_root):
    """Write batch to Parquet, grouped by dataset_id."""
    if not rows:
        return

    df = pd.DataFrame(rows)

    # Group by dataset_id to avoid partition limit
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


def split_fasta_into_chunks(fasta_path, num_chunks):
    """
    Read FASTA and split into chunks.

    Args:
        fasta_path: Path to FASTA file
        num_chunks: Number of chunks to create

    Returns:
        List of record lists (chunks)
    """
    print(f"\n📖 Reading FASTA file: {fasta_path}...")
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total = len(records)
    print(f"✓ Loaded {total:,} sequences")

    # Split into chunks
    chunk_size = (total + num_chunks - 1) // num_chunks  # Ceiling division
    chunks = []

    for i in range(0, total, chunk_size):
        chunk = records[i : i + chunk_size]
        chunks.append(chunk)

    print(f"✓ Split into {len(chunks)} chunks (~{chunk_size:,} sequences each)")
    return chunks


def parallel_fasta_to_parquet(
    fasta_path: Path,
    metadata_dict: dict[str, dict[str, object]],
    skani_rep_dict: dict[str, bool],
    parquet_root: Path,
    orthogroup_dict: dict[str, str] | None = None,
    singleton_set: set[str] | None = None,
    eggnog_dict: dict[str, dict[str, object | None]] | None = None,
    num_workers: int = 16,
) -> None:
    """
    Ingest FASTA in parallel using multiprocessing.

    Args:
        fasta_path: Path to FASTA file
        metadata_dict: Genome metadata dictionary
        skani_rep_dict: Skani representative status dictionary
        parquet_root: Output directory
        orthogroup_dict: Protein → orthogroup_id mapping (optional)
        singleton_set: Set of singleton protein IDs (optional)
        eggnog_dict: Protein → eggNOG annotations mapping (optional)
        num_workers: Number of parallel workers
    """
    start_time = time.time()

    # Split FASTA into chunks
    chunks = split_fasta_into_chunks(fasta_path, num_workers)

    # Prepare arguments for each worker
    # Annotations are passed to all workers via copy-on-write (Linux)
    args_list = [
        (
            chunk,
            metadata_dict,
            skani_rep_dict,
            parquet_root,
            orthogroup_dict,
            singleton_set,
            eggnog_dict,
            i + 1,
        )
        for i, chunk in enumerate(chunks)
    ]

    # Process in parallel
    print(f"\n🚀 Processing with {num_workers} workers...")
    with Pool(processes=num_workers) as pool:
        results = pool.map(process_fasta_chunk, args_list)

    # Aggregate results
    total_processed = sum(r[0] for r in results)
    total_skipped = sum(r[1] for r in results)

    elapsed = time.time() - start_time

    print("\n✓ Parallel ingestion complete:")
    print(f"  - Processed: {total_processed:,} sequences")
    print(f"  - Skipped (no metadata): {total_skipped:,}")
    print(f"  - Time: {elapsed/60:.1f} minutes")
    print(f"  - Speed: {total_processed/elapsed:.0f} sequences/second")
    print(f"  - Output: {parquet_root}")


def ingest_all_parallel(
    metadata_path: Path,
    fna_path: Path,
    faa_path: Path,
    parquet_root: Path = Path("data/parquet/sequences"),
    num_workers: int = 16,
):
    """
    Ingest both nucleotide and protein sequences in parallel with annotations.

    Args:
        metadata_path: Path to metadata TSV
        fna_path: Path to nucleotide FASTA
        faa_path: Path to protein FASTA
        parquet_root: Root for Parquet output
        num_workers: Number of parallel workers (default: 16)
    """
    parquet_root.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("GVMAG Database Parallel Ingestion")
    print("=" * 70)
    print(f"Workers: {num_workers} CPUs")

    # Load metadata BEFORE forking workers (for copy-on-write)
    print(f"\n📖 Loading metadata from {metadata_path}...")
    metadata_dict, skani_rep_dict = parse_metadata_tsv(metadata_path)
    print(f"✓ Loaded metadata for {len(metadata_dict):,} genomes")
    print(f"  - Skani representatives: {sum(skani_rep_dict.values())}")

    # Load annotations BEFORE forking workers (for copy-on-write memory sharing)
    print("\n📊 Loading protein annotations...")
    # Use project_root defined at module level
    root = Path(__file__).parent.parent
    proteinortho_tsv = root / "ingestion_data/proteinortho/gvmags_db_trial4.proteinortho.tsv"
    singletons_tsv = root / "ingestion_data/proteinortho/gvmags_db_trial4.singletons.tsv"
    eggnog_tsv = root / "ingestion_data/annotation/emapper_out.emapper.annotations"

    orthogroup_dict = None
    singleton_set = None
    eggnog_dict = None

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
    parallel_fasta_to_parquet(
        fna_path,
        metadata_dict,
        skani_rep_dict,
        parquet_root,
        orthogroup_dict,
        singleton_set,
        eggnog_dict,
        num_workers,
    )

    # Ingest proteins
    print("\n🧬 Ingesting protein sequences...")
    parallel_fasta_to_parquet(
        faa_path,
        metadata_dict,
        skani_rep_dict,
        parquet_root,
        orthogroup_dict,
        singleton_set,
        eggnog_dict,
        num_workers,
    )

    print("\n" + "=" * 70)
    print("✓ All ingestion complete!")
    print("=" * 70)


if __name__ == "__main__":
    # Configuration - use paths relative to project root
    root = Path(__file__).parent.parent
    metadata = root / "ingestion_data/Updated_naming_Sept2025.tsv"
    fna = root / "ingestion_data/gvmagsV2all.fna"
    faa = root / "ingestion_data/gvmagsV2all.faa"
    output = root / "artifacts/parquet/sequences"
    workers = 16  # Use 16 CPUs

    try:
        ingest_all_parallel(metadata, fna, faa, output, num_workers=workers)
    except Exception as e:
        print(f"\n❌ Ingestion FAILED: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
