"""
Comprehensive genome placement analysis.

Combines TNF distance and skani ANI to determine taxonomic placement
of query genomes against the database.
"""

import json
from pathlib import Path

from gvmagdb.analyze.run_skani import skani_with_metadata
from gvmagdb.analyze.tnf_distance import find_closest_by_tnf


def analyze_genome(
    query_fna: Path,
    ref_fna: Path | None = None,
    output_dir: Path | None = None,
    distance_metric: str = "cosine",
    top_k: int = 20,
    ani_threshold: float = 80.0,
    threads: int = 8,
    skip_skani: bool = False,
    skip_tnf: bool = False,
) -> dict:
    """
    Comprehensive genome placement analysis using TNF and ANI.

    Args:
        query_fna: Query genome FASTA
        ref_fna: Reference genome database (default: artifacts/products/fna/gvmagdb.fna)
        output_dir: Output directory (default: {query_basename}_analysis/)
        distance_metric: TNF distance metric (cosine or euclidean)
        top_k: Number of top TNF matches to return
        ani_threshold: Minimum ANI to report (default: 80%)
        threads: Number of threads
        skip_skani: Skip skani ANI calculation
        skip_tnf: Skip TNF distance calculation

    Returns:
        Dictionary with placement results
    """
    print("\n" + "=" * 70)
    print("GENOME PLACEMENT ANALYSIS")
    print("=" * 70)
    print(f"Query: {query_fna}")

    # Set defaults
    if ref_fna is None:
        ref_fna = Path("artifacts/products/fna/gvmagdb.fna")

    if output_dir is None:
        query_basename = query_fna.stem
        output_dir = Path(f"{query_basename}_analysis")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "query": str(query_fna),
        "tnf_results": None,
        "ani_results": None,
        "summary": {},
    }

    # TNF distance calculation
    if not skip_tnf:
        print("\n" + "=" * 70)
        print("TNF DISTANCE ANALYSIS")
        print("=" * 70)

        tnf_results = find_closest_by_tnf(
            query_fasta=query_fna,
            distance_metric=distance_metric,
            top_k=top_k,
        )

        # Save TNF results
        tnf_output = output_dir / "tnf_distances.tsv"
        tnf_results.to_csv(tnf_output, sep="\t", index=False)
        print(f"\n✓ TNF results saved: {tnf_output}")

        results["tnf_results"] = tnf_results
        results["summary"]["tnf_top_match"] = {
            "dataset_id": tnf_results.iloc[0]["dataset_id"],
            "distance": float(tnf_results.iloc[0]["tnf_distance"]),
            "taxonomy": tnf_results.iloc[0].get("taxonomy_majority", "N/A"),
        }

        # Show top 5 matches
        print("\n🔝 Top 5 TNF matches:\n")
        display_cols = ["dataset_id", "tnf_distance", "genome_length", "phylum", "s_cluster"]
        available_cols = [col for col in display_cols if col in tnf_results.columns]
        print(tnf_results.head(5)[available_cols].to_string(index=False))

    # Skani ANI calculation
    if not skip_skani:
        if ref_fna.exists():
            print("\n" + "=" * 70)
            print("SKANI ANI ANALYSIS")
            print("=" * 70)

            ani_results = skani_with_metadata(
                query_fna, ref_fna, ani_threshold=ani_threshold, threads=threads
            )

            if not ani_results.empty:
                # Save ANI results
                ani_output = output_dir / "ani_results.tsv"
                ani_results.to_csv(ani_output, sep="\t", index=False)
                print(f"\n✓ ANI results saved: {ani_output}")

                results["ani_results"] = ani_results
                results["summary"]["ani_top_match"] = {
                    "dataset_id": ani_results.iloc[0]["Ref_name"],
                    "ani": float(ani_results.iloc[0]["ANI"]),
                    "cluster": ani_results.iloc[0].get("s_cluster", "N/A"),
                }

                # Show top 5 matches
                print("\n🔝 Top 5 ANI matches:\n")
                display_cols = ["Ref_name", "ANI", "Align_fraction_query", "s_cluster", "genome"]
                available_cols = [col for col in display_cols if col in ani_results.columns]
                print(ani_results.head(5)[available_cols].to_string(index=False))
            else:
                print(f"\n⚠️  No ANI hits found above {ani_threshold}% threshold")
                results["summary"]["ani_message"] = f"No hits ≥{ani_threshold}% ANI"
        else:
            print(f"\n⚠️  Reference database not found: {ref_fna}")
            print(
                "   Export it with: gvmagdb export-fna-cmd --output artifacts/products/fna/gvmagdb.fna"
            )
            results["summary"]["ani_message"] = "Reference database not found"

    # Generate summary report
    summary_output = output_dir / "placement_summary.txt"
    _write_summary_report(results, summary_output)
    print(f"\n✓ Summary report: {summary_output}")

    # Save JSON results
    json_output = output_dir / "placement_results.json"
    _save_json_results(results, json_output)
    print(f"✓ JSON results: {json_output}")

    print(f"\n✓ Analysis complete! Output directory: {output_dir}/")

    return results


def _write_summary_report(results: dict, output_path: Path):
    """Write human-readable summary report."""
    with open(output_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("GENOME PLACEMENT ANALYSIS SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Query: {results['query']}\n\n")

        # TNF results
        if results.get("tnf_results") is not None:
            f.write("TNF DISTANCE RESULTS\n")
            f.write("-" * 70 + "\n")
            tnf_top = results["summary"]["tnf_top_match"]
            f.write(f"Top match: {tnf_top['dataset_id']}\n")
            f.write(f"Distance: {tnf_top['distance']:.6f}\n")
            f.write(f"Taxonomy: {tnf_top['taxonomy']}\n\n")

        # ANI results
        if results.get("ani_results") is not None:
            f.write("SKANI ANI RESULTS\n")
            f.write("-" * 70 + "\n")
            ani_top = results["summary"]["ani_top_match"]
            f.write(f"Top match: {ani_top['dataset_id']}\n")
            f.write(f"ANI: {ani_top['ani']:.2f}%\n")
            f.write(f"Cluster: {ani_top['cluster']}\n\n")
        elif "ani_message" in results["summary"]:
            f.write("SKANI ANI RESULTS\n")
            f.write("-" * 70 + "\n")
            f.write(f"{results['summary']['ani_message']}\n\n")

        # Interpretation
        f.write("INTERPRETATION\n")
        f.write("-" * 70 + "\n")

        if results.get("ani_results") is not None and len(results["ani_results"]) > 0:
            top_ani = results["summary"]["ani_top_match"]["ani"]
            if top_ani >= 95:
                f.write("Species-level match found (≥95% ANI)\n")
                f.write(
                    f"Query belongs to cluster: {results['summary']['ani_top_match']['cluster']}\n"
                )
            elif top_ani >= 90:
                f.write("Genus-level similarity (90-95% ANI)\n")
                f.write("Query likely represents a related but distinct species\n")
            else:
                f.write("Distant relationship (<90% ANI)\n")
                f.write("Query may represent a novel species or genus\n")
        else:
            f.write("No close ANI matches found - likely novel lineage\n")

        if results.get("tnf_results") is not None:
            tnf_dist = results["summary"]["tnf_top_match"]["distance"]
            f.write(f"\nTNF compositional similarity: {1-tnf_dist:.3f}\n")
            taxonomy_str = results["summary"]["tnf_top_match"]["taxonomy"]
            if taxonomy_str != "N/A":
                f.write(f"Closest compositional match: {taxonomy_str}\n")


def _save_json_results(results: dict, output_path: Path):
    """Save machine-readable JSON results."""
    # Convert DataFrames to dicts for JSON serialization
    json_data = {
        "query": results["query"],
        "summary": results["summary"],
    }

    if results.get("tnf_results") is not None:
        json_data["tnf_top_matches"] = results["tnf_results"].head(10).to_dict(orient="records")

    if results.get("ani_results") is not None:
        json_data["ani_top_matches"] = results["ani_results"].head(10).to_dict(orient="records")

    with open(output_path, "w") as f:
        json.dump(json_data, f, indent=2, default=str)
