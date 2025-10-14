"""
Skani ANI search with join-back to Parquet.

Searches query genome against reference database and joins results
back to genome metadata including cluster information.
"""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from gvmagdb.core.catalog import connect

# Skani search output columns
SKANI_COLUMNS = [
    "Ref_file",
    "Query_file",
    "ANI",
    "Align_fraction_ref",
    "Align_fraction_query",
    "Ref_name",
    "Query_name",
]


def run_skani_search(
    query_fna: Path,
    ref_fna: Path,
    output_tsv: Path,
    ani_threshold: float = 80.0,
    marker_compression: int = 1000,
    threads: int = 8,
) -> pd.DataFrame:
    """
    Run skani search for ANI calculation.

    Args:
        query_fna: Query genome FASTA
        ref_fna: Reference database FASTA
        output_tsv: Output TSV path
        ani_threshold: Minimum ANI to report (default 80%)
        marker_compression: Marker k-mer compression factor (default 1000, skani default)
        threads: Number of threads

    Returns:
        DataFrame with skani search results
    """
    print("\n🔍 Running skani ANI calculation...")
    print(f"   Query: {query_fna}")
    print(f"   Reference: {ref_fna}")
    print(f"   ANI threshold: {ani_threshold}%")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Skani dist command for FASTA-to-FASTA comparison
    # Using dist mode for direct comparison without pre-built database
    cmd = [
        "skani",
        "dist",
        str(query_fna),
        str(ref_fna),
        "-o",
        str(output_tsv),
        "-t",
        str(threads),
        "-m",
        str(marker_compression),  # Marker compression for filtering
        "--ci",  # Confidence intervals
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # Print skani info messages (they go to stderr)
        if result.stderr:
            for line in result.stderr.split("\n"):
                if "INFO" in line or "WARN" in line:
                    print(f"   {line.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Skani failed: {e.stderr}")
        raise

    # Parse results
    if output_tsv.exists() and output_tsv.stat().st_size > 0:
        df = pd.read_csv(output_tsv, sep="\t")

        # Filter by ANI threshold
        df = df[df["ANI"] >= ani_threshold].copy()

        print(f"✓ Found {len(df)} hits above {ani_threshold}% ANI")
        return df
    else:
        print("✓ No hits found")
        return pd.DataFrame(columns=SKANI_COLUMNS)


def skani_with_metadata(
    query_fna: Path,
    ref_fna: Path,
    ani_threshold: float = 80.0,
    marker_compression: int = 1000,
    threads: int = 8,
) -> pd.DataFrame:
    """
    Run skani search and join with genome metadata including cluster information.

    Args:
        query_fna: Query genome FASTA
        ref_fna: Reference database FASTA
        ani_threshold: Minimum ANI to report (default 80%)
        marker_compression: Marker k-mer compression factor (default 1000, skani default)
        threads: Number of threads

    Returns:
        DataFrame with skani results + genome metadata + cluster info
    """
    # Run skani search
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
        output_tsv = Path(tmp.name)

    try:
        results = run_skani_search(
            query_fna, ref_fna, output_tsv, ani_threshold, marker_compression, threads
        )

        if results.empty:
            return results

        print("\n🔗 Joining results with genome metadata and clusters...")

        # Connect to catalog
        conn = connect()

        try:
            # Get unique reference genome IDs (Ref_name should match dataset_id)
            ref_genomes = results["Ref_name"].unique().tolist()

            # Build placeholders for SQL IN clause
            placeholders = ",".join(["?"] * len(ref_genomes))

            # Query genome metadata including cluster information
            # Group by dataset_id since we want genome-level info
            query = f"""
                SELECT DISTINCT
                    dataset_id,
                    common_name as genome,
                    phylum,
                    class,
                    order,
                    family,
                    genus,
                    gvclass_phylum,
                    gvclass_class,
                    gvclass_order,
                    gvclass_family,
                    gvclass_genus,
                    gvclass_species,
                    taxonomy_majority,
                    s_cluster,
                    skani_rep,
                    pdm_cluster,
                    mcl_rep,
                    lenbp,
                    gcperc,
                    genecount,
                    order_completeness,
                    order_confidence_score,
                    ecosystem,
                    habitat
                FROM sequences
                WHERE dataset_id IN ({placeholders})
            """

            metadata = conn.execute(query, ref_genomes).df()

            # Join results with metadata
            annotated = results.merge(
                metadata, left_on="Ref_name", right_on="dataset_id", how="left"
            )

            # Sort by ANI (highest first)
            annotated = annotated.sort_values("ANI", ascending=False)

            print(f"✓ Annotated {len(annotated)} hits with metadata")

            # Show cluster summary
            if "s_cluster" in annotated.columns:
                clusters = annotated["s_cluster"].value_counts()
                print("\n📊 Cluster distribution:")
                for cluster, count in clusters.head(5).items():
                    print(f"   {cluster}: {count} genomes")

            return annotated

        finally:
            conn.close()

    finally:
        # Clean up temp file
        if output_tsv.exists():
            output_tsv.unlink()


def summarize_cluster_membership(results: pd.DataFrame, species_threshold: float = 95.0) -> dict:
    """
    Summarize which clusters the query genome belongs to.

    Args:
        results: Annotated skani results DataFrame
        species_threshold: ANI threshold for species-level clustering (default 95%)

    Returns:
        Dictionary with cluster membership summary
    """
    if results.empty:
        return {"message": "No matches found"}

    # Filter for species-level matches (≥95% ANI)
    species_matches = results[results["ANI"] >= species_threshold].copy()

    summary = {
        "total_hits": len(results),
        "species_level_matches": len(species_matches),
        "top_match_ani": float(results.iloc[0]["ANI"]) if len(results) > 0 else 0.0,
        "top_match_genome": results.iloc[0]["genome"] if len(results) > 0 else "N/A",
    }

    # Cluster membership
    if len(species_matches) > 0:
        # Most common s_cluster
        if "s_cluster" in species_matches.columns:
            cluster_counts = species_matches["s_cluster"].value_counts()
            summary["likely_species_cluster"] = (
                cluster_counts.index[0] if len(cluster_counts) > 0 else "Unknown"
            )
            summary["species_cluster_genomes"] = (
                int(cluster_counts.iloc[0]) if len(cluster_counts) > 0 else 0
            )

        # Representative genome
        if "skani_rep" in species_matches.columns:
            summary["cluster_representative"] = species_matches.iloc[0]["skani_rep"]

        # Phylogenetic cluster
        if "pdm_cluster" in species_matches.columns:
            pdm_counts = species_matches["pdm_cluster"].value_counts()
            summary["phylogenetic_cluster"] = (
                pdm_counts.index[0] if len(pdm_counts) > 0 else "Unknown"
            )
    else:
        summary["message"] = (
            f"No species-level matches (≥{species_threshold}% ANI). Top hit: {summary['top_match_ani']:.2f}%"
        )

    return summary
