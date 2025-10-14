"""
PyHMMER search with join-back to Parquet.

Searches HMM profile against protein sequences and joins results
back to sequence metadata via gid.
"""

import tempfile
from pathlib import Path

import pandas as pd
import pyhmmer

from gvmagdb.core.catalog import connect
from gvmagdb.export.export_faa import export_aa_fasta


def hmmsearch(
    hmm_path: Path,
    target_faa: Path,
    output_tsv: Path,
    evalue: float = 1e-5,
    threads: int = 8,
) -> pd.DataFrame:
    """
    Run PyHMMER hmmsearch.

    Args:
        hmm_path: HMM profile file
        target_faa: Target protein FASTA
        output_tsv: Output TSV path
        evalue: E-value threshold
        threads: Number of threads

    Returns:
        DataFrame with HMM search results
    """
    print("\n🔍 Running PyHMMER search...")
    print(f"   HMM: {hmm_path}")
    print(f"   Target: {target_faa}")
    print(f"   E-value: {evalue}")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Load HMM profile
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        hmm = hmm_file.read()

    # Load target sequences
    with pyhmmer.easel.SequenceFile(
        target_faa, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()
    ) as seq_file:
        sequences = list(seq_file)

    print(f"   Loaded {len(sequences)} target sequences")

    # Run search
    results = []
    for hits in pyhmmer.hmmsearch([hmm], sequences, cpus=threads):
        for hit in hits:
            if hit.evalue <= evalue:
                results.append(
                    {
                        "query_name": (
                            hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                        ),
                        "target_name": (
                            hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                        ),
                        "evalue": hit.evalue,
                        "score": hit.score,
                        "bias": hit.bias,
                    }
                )

    if results:
        df = pd.DataFrame(results)
        df.to_csv(output_tsv, sep="\t", index=False)
        print(f"✓ Found {len(df)} hits")
        return df
    else:
        print("✓ No hits found")
        return pd.DataFrame(columns=["query_name", "target_name", "evalue", "score", "bias"])


def hmmsearch_with_metadata(
    hmm_path: Path,
    target_faa: Path | None = None,
    dataset_id: str | None = None,
    evalue: float = 1e-5,
    threads: int = 8,
) -> pd.DataFrame:
    """
    Run PyHMMER hmmsearch and join with sequence metadata.

    Args:
        hmm_path: HMM profile file
        target_faa: Target protein FASTA (if None, exports from Parquet)
        dataset_id: Filter by dataset_id (None = all)
        evalue: E-value threshold
        threads: Number of threads

    Returns:
        DataFrame with HMM results + metadata
    """
    # Export target sequences if not provided
    if target_faa is None:
        print("\n📤 Exporting target sequences...")
        temp_faa = Path(tempfile.mktemp(suffix=".faa"))
        export_aa_fasta(temp_faa, dataset_id=dataset_id)
        target_faa = temp_faa
        cleanup = True
    else:
        cleanup = False

    try:
        # Run search
        output_tsv = Path("results/hmm_results.tsv")
        results = hmmsearch(hmm_path, target_faa, output_tsv, evalue, threads)

        if results.empty:
            return results

        print("\n🔗 Joining results with metadata...")

        # Connect to catalog
        conn = connect()

        try:
            # Get unique target gids (target_name IS gid from our export)
            target_gids = results["target_name"].unique().tolist()

            # Build placeholders for SQL IN clause
            placeholders = ",".join(["?"] * len(target_gids))

            # Query metadata for all target sequences
            query = f"""
                SELECT
                    gid,
                    dataset_id,
                    length,
                    genome,
                    phylum,
                    class,
                    order,
                    family,
                    genus,
                    ecosystem,
                    habitat
                FROM sequences
                WHERE gid IN ({placeholders}) AND seq_type = 'AA'
            """

            metadata = conn.execute(query, target_gids).df()

            # Join results with metadata
            annotated = results.merge(metadata, left_on="target_name", right_on="gid", how="left")

            print(f"✓ Annotated {len(annotated)} hits with metadata")
            return annotated

        finally:
            conn.close()

    finally:
        if cleanup and target_faa.exists():
            target_faa.unlink()
