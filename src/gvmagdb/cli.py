"""
Command-line interface for GVMAG Database.
"""

import sys
from pathlib import Path

import click

from gvmagdb.analyze import analyze_genome, run_blastp, run_pyhmmer, run_skani
from gvmagdb.core import catalog
from gvmagdb.export import build_dmnd, export_faa, export_fna, extract_hits
from gvmagdb.prepare import ingest


@click.group()
def cli():
    """GVMAG Genome Database - instructions.md compliant implementation."""
    pass


@cli.command()
@click.option(
    "--fna",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Nucleotide FASTA file",
)
@click.option(
    "--faa", type=click.Path(exists=True, path_type=Path), required=True, help="Protein FASTA file"
)
@click.option(
    "--metadata",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Metadata TSV file",
)
@click.option(
    "--output",
    type=click.Path(path_type=Path),
    default=Path("artifacts/parquet/sequences"),
    help="Parquet output directory",
)
def ingest_cmd(fna, faa, metadata, output):
    """Ingest FASTA files to Parquet."""
    try:
        ingest.ingest_all(metadata, fna, faa, output)
    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--parquet",
    type=click.Path(path_type=Path),
    default=Path("artifacts/parquet/sequences/**/*.parquet"),
    help="Parquet glob pattern",
)
def stats(parquet):
    """Show database statistics."""
    try:
        conn = catalog.connect(parquet_glob=str(parquet))
        stats_data = catalog.get_stats(conn)

        click.echo("\n" + "=" * 70)
        click.echo("Database Statistics")
        click.echo("=" * 70)

        click.echo(f"\n📊 Total sequences: {stats_data['total_sequences']:,}")
        click.echo(f"   Unique genomes: {stats_data['unique_genomes']:,}")
        click.echo(f"   Total bases: {stats_data['total_bases']:,}")

        click.echo("\n🧬 By sequence type:")
        for seq_type, count in stats_data["by_seq_type"].items():
            click.echo(f"   {seq_type}: {count:,}")

        click.echo(f"\n📈 Average GC content (NT): {stats_data['avg_gc_content']:.2f}%")

        click.echo("\n🔝 Top 10 datasets:")
        for dataset_id, count in list(stats_data["top_10_datasets"].items())[:10]:
            click.echo(f"   {dataset_id}: {count:,}")

        click.echo("\n" + "=" * 70 + "\n")

        conn.close()

    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--output",
    type=click.Path(path_type=Path),
    default=Path("artifacts/products/faa/gvmagdb.faa"),
    help="Output FASTA file (default: artifacts/products/faa/gvmagdb.faa)",
)
@click.option("--dataset", type=str, default=None, help="Filter by dataset_id")
def export_faa_cmd(output, dataset):
    """Export protein sequences to FASTA."""
    try:
        export_faa.export_aa_fasta(output, dataset_id=dataset)
    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--output", type=click.Path(path_type=Path), required=True, help="Output FASTA file")
@click.option("--dataset", type=str, default=None, help="Filter by dataset_id")
def export_fna_cmd(output, dataset):
    """Export nucleotide sequences to FASTA."""
    try:
        export_fna.export_nt_fasta(output, dataset_id=dataset)
    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--output",
    type=click.Path(path_type=Path),
    default=Path("artifacts/products/dmnd/gvmagdb.dmnd"),
    help="Output DIAMOND database (default: artifacts/products/dmnd/gvmagdb.dmnd)",
)
@click.option("--dataset", type=str, default=None, help="Filter by dataset_id (None = all)")
@click.option("--threads", type=int, default=8, help="Number of threads")
def build_diamond(output, dataset, threads):
    """Build DIAMOND database from gvmagDB sequences."""
    try:
        build_dmnd.build_diamond_db(dataset, output, threads)
    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--query",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Query protein FASTA",
)
@click.option(
    "--db",
    type=click.Path(path_type=Path),
    default=None,
    help="DIAMOND database (default: artifacts/products/dmnd/gvmagdb.dmnd)",
)
@click.option(
    "--outdir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory (default: {query_basename}_vs_gvmagdb/)",
)
@click.option("--evalue", type=float, default=1e-5, help="E-value threshold")
@click.option("--threads", type=int, default=8, help="Number of threads")
@click.option("--extract/--no-extract", default=True, help="Extract hit sequences (default: True)")
def diamond(query, db, outdir, evalue, threads, extract):
    """Run DIAMOND BLASTP search with metadata join and optional hit extraction."""
    try:
        # Default database path
        if db is None:
            db = Path("artifacts/products/dmnd/gvmagdb.dmnd")
            if not db.exists():
                click.echo(f"\n❌ Default database not found: {db}", err=True)
                click.echo("   Build it with: gvmagdb build-diamond", err=True)
                sys.exit(1)
        else:
            db = Path(db)
            if not db.exists():
                click.echo(f"\n❌ Database not found: {db}", err=True)
                sys.exit(1)

        # Default output directory
        if outdir is None:
            query_basename = Path(query).stem
            outdir = Path(f"{query_basename}_vs_gvmagdb")

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Output files
        m8_output = outdir / "gvmagdb_hits.m8"
        tsv_output = outdir / "gvmagdb_hits.tsv"
        faa_output = outdir / "gvmagdb_hits.faa"

        # Run search
        results = run_blastp.blastp_with_metadata(query, db, evalue=evalue, threads=threads)

        if not results.empty:
            # Save BLAST tabular format (m8)
            blast_cols = [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ]
            results[blast_cols].to_csv(m8_output, sep="\t", index=False, header=False)
            click.echo(f"\n✓ BLAST results: {m8_output}")

            # Save full results with metadata
            results.to_csv(tsv_output, sep="\t", index=False)
            click.echo(f"✓ Annotated results: {tsv_output}")

            # Extract hit sequences
            if extract:
                extract_hits.extract_sequences(
                    hits_tsv=tsv_output, output_fasta=faa_output, gid_column="sseqid", seq_type="AA"
                )
                click.echo(f"✓ Extracted sequences: {faa_output}")

            # Show summary
            click.echo("\n📊 Summary:")
            click.echo(f"   Total hits: {len(results)}")
            click.echo(f"   Unique queries with hits: {results['qseqid'].nunique()}")
            click.echo(f"   Output directory: {outdir}/")

            # Show top 10 hits
            click.echo("\n🔝 Top 10 hits:\n")
            top10 = results.head(10)[
                ["qseqid", "sseqid", "pident", "evalue", "bitscore", "phylum", "genome"]
            ]
            click.echo(top10.to_string(index=False))
            click.echo()

        else:
            click.echo("\n⚠ No hits found")

    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--hmm", type=click.Path(exists=True, path_type=Path), required=True, help="HMM profile file"
)
@click.option(
    "--target",
    type=click.Path(path_type=Path),
    default=None,
    help="Target protein FASTA (default: artifacts/products/faa/gvmagdb.faa)",
)
@click.option(
    "--outdir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory (default: {hmm_basename}_vs_gvmagdb/)",
)
@click.option("--evalue", type=float, default=1e-5, help="E-value threshold")
@click.option("--threads", type=int, default=8, help="Number of threads")
@click.option("--extract/--no-extract", default=True, help="Extract hit sequences (default: True)")
def hmmsearch(hmm, target, outdir, evalue, threads, extract):
    """Run PyHMMER search with metadata join and optional hit extraction."""
    try:
        # Default target database
        if target is None:
            target = Path("artifacts/products/faa/gvmagdb.faa")
            if not target.exists():
                click.echo(f"\n❌ Default target not found: {target}", err=True)
                click.echo(
                    "   Build it with: gvmagdb export-faa-cmd --output artifacts/products/faa/gvmagdb.faa",
                    err=True,
                )
                sys.exit(1)
        else:
            target = Path(target)
            if not target.exists():
                click.echo(f"\n❌ Target not found: {target}", err=True)
                sys.exit(1)

        # Default output directory
        if outdir is None:
            hmm_basename = Path(hmm).stem
            outdir = Path(f"{hmm_basename}_vs_gvmagdb")

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Output files
        tbl_output = outdir / "gvmagdb_hits.tbl"
        tsv_output = outdir / "gvmagdb_hits.tsv"
        faa_output = outdir / "gvmagdb_hits.faa"

        # Run search
        results = run_pyhmmer.hmmsearch_with_metadata(
            hmm, target_faa=target, dataset_id=None, evalue=evalue, threads=threads
        )

        if not results.empty:
            # Save HMM tabular format (simplified)
            hmm_cols = ["query_name", "target_name", "evalue", "score", "bias"]
            results[hmm_cols].to_csv(tbl_output, sep="\t", index=False, header=False)
            click.echo(f"\n✓ HMM results: {tbl_output}")

            # Save full results with metadata
            results.to_csv(tsv_output, sep="\t", index=False)
            click.echo(f"✓ Annotated results: {tsv_output}")

            # Extract hit sequences
            if extract:
                extract_hits.extract_sequences(
                    hits_tsv=tsv_output,
                    output_fasta=faa_output,
                    gid_column="target_name",
                    seq_type="AA",
                )
                click.echo(f"✓ Extracted sequences: {faa_output}")

            # Show summary
            click.echo("\n📊 Summary:")
            click.echo(f"   Total hits: {len(results)}")
            click.echo(f"   Unique HMMs with hits: {results['query_name'].nunique()}")
            click.echo(f"   Output directory: {outdir}/")

            # Show top 10 hits
            click.echo("\n🔝 Top 10 hits:\n")
            top10 = results.head(10)[
                ["query_name", "target_name", "evalue", "score", "phylum", "genome"]
            ]
            click.echo(top10.to_string(index=False))
            click.echo()

        else:
            click.echo("\n⚠ No hits found")

    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--query",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Query genome FASTA (nucleotides)",
)
@click.option(
    "--ref",
    type=click.Path(path_type=Path),
    default=None,
    help="Reference genome database (default: artifacts/products/fna/gvmagdb.fna)",
)
@click.option(
    "--outdir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory (default: {query_basename}_vs_gvmagdb/)",
)
@click.option(
    "--ani-threshold", type=float, default=80.0, help="Minimum ANI to report (default: 80%)"
)
@click.option(
    "--species-threshold",
    type=float,
    default=95.0,
    help="ANI threshold for species-level clustering (default: 95%)",
)
@click.option(
    "--marker-compression",
    type=int,
    default=1000,
    help="Marker k-mer compression factor (default: 1000, skani default)",
)
@click.option("--threads", type=int, default=8, help="Number of threads")
def skani(query, ref, outdir, ani_threshold, species_threshold, marker_compression, threads):
    """Run skani ANI search to find genome cluster membership."""
    try:
        # Default reference database
        if ref is None:
            ref = Path("artifacts/products/fna/gvmagdb.fna")
            if not ref.exists():
                click.echo(f"\n❌ Default reference not found: {ref}", err=True)
                click.echo(
                    "   Build it with: gvmagdb export-fna-cmd --output artifacts/products/fna/gvmagdb.fna",
                    err=True,
                )
                sys.exit(1)
        else:
            ref = Path(ref)
            if not ref.exists():
                click.echo(f"\n❌ Reference not found: {ref}", err=True)
                sys.exit(1)

        # Default output directory
        if outdir is None:
            query_basename = Path(query).stem
            outdir = Path(f"{query_basename}_vs_gvmagdb")

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Output files
        tsv_output = outdir / "gvmagdb_ani.tsv"
        summary_output = outdir / "cluster_membership.txt"

        # Run skani search
        results = run_skani.skani_with_metadata(
            query,
            ref,
            ani_threshold=ani_threshold,
            marker_compression=marker_compression,
            threads=threads,
        )

        if not results.empty:
            # Save full results with metadata
            results.to_csv(tsv_output, sep="\t", index=False)
            click.echo(f"\n✓ ANI results: {tsv_output}")

            # Summarize cluster membership
            summary = run_skani.summarize_cluster_membership(results, species_threshold)

            # Save summary
            with open(summary_output, "w") as f:
                f.write("=" * 70 + "\n")
                f.write("Cluster Membership Summary\n")
                f.write("=" * 70 + "\n\n")
                for key, value in summary.items():
                    f.write(f"{key}: {value}\n")

            click.echo(f"✓ Cluster summary: {summary_output}")

            # Show summary on screen
            click.echo("\n📊 Cluster Membership Summary:")
            click.echo(f"   Total hits: {summary.get('total_hits', 0)}")
            click.echo(
                f"   Species-level matches (≥{species_threshold}% ANI): {summary.get('species_level_matches', 0)}"
            )
            click.echo(
                f"   Top match: {summary.get('top_match_genome', 'N/A')} ({summary.get('top_match_ani', 0):.2f}% ANI)"
            )

            if "likely_species_cluster" in summary:
                click.echo(f"\n🎯 Likely Species Cluster: {summary['likely_species_cluster']}")
                click.echo(f"   Cluster has {summary.get('species_cluster_genomes', 0)} genome(s)")
                if "cluster_representative" in summary:
                    click.echo(f"   Representative: {summary['cluster_representative']}")
                if "phylogenetic_cluster" in summary:
                    click.echo(f"   Phylogenetic cluster: {summary['phylogenetic_cluster']}")
            elif "message" in summary:
                click.echo(f"\n⚠️  {summary['message']}")

            click.echo(f"\n   Output directory: {outdir}/")

            # Show top 10 hits
            click.echo("\n🔝 Top 10 ANI hits:\n")
            display_cols = [
                "Query_name",
                "Ref_name",
                "ANI",
                "Align_fraction_query",
                "s_cluster",
                "genome",
            ]
            available_cols = [col for col in display_cols if col in results.columns]
            top10 = results.head(10)[available_cols]
            click.echo(top10.to_string(index=False))
            click.echo()

        else:
            click.echo(f"\n⚠ No hits found above {ani_threshold}% ANI threshold")

    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        import traceback

        traceback.print_exc()
        sys.exit(1)


@cli.command()
@click.option(
    "--hits",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Search results TSV with gid column",
)
@click.option("--fasta", type=click.Path(path_type=Path), default=None, help="Output FASTA file")
@click.option(
    "--tsv",
    type=click.Path(path_type=Path),
    default=None,
    help="Output TSV file with full metadata",
)
@click.option(
    "--json",
    type=click.Path(path_type=Path),
    default=None,
    help="Output JSON file with full metadata",
)
@click.option(
    "--gid-column", type=str, default="sseqid", help="Column name containing gids (default: sseqid)"
)
@click.option(
    "--seq-type", type=click.Choice(["AA", "NT"]), default=None, help="Filter by sequence type"
)
def extract(hits, fasta, tsv, json, gid_column, seq_type):
    """Extract sequences from search results."""
    try:
        if not any([fasta, tsv, json]):
            click.echo(
                "\n❌ Error: Specify at least one output format (--fasta, --tsv, or --json)",
                err=True,
            )
            sys.exit(1)

        extract_hits.extract_sequences(
            hits,
            output_fasta=fasta,
            output_tsv=tsv,
            output_json=json,
            gid_column=gid_column,
            seq_type=seq_type,
        )

        click.echo("\n✓ Extraction complete")

    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--query",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Query genome FASTA (nucleotides)",
)
@click.option(
    "--ref",
    type=click.Path(path_type=Path),
    default=None,
    help="Reference genome database (default: artifacts/products/fna/gvmagdb.fna)",
)
@click.option(
    "--outdir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory (default: {query_basename}_analysis/)",
)
@click.option(
    "--distance-metric",
    type=click.Choice(["cosine", "euclidean"]),
    default="cosine",
    help="TNF distance metric (default: cosine)",
)
@click.option(
    "--top-k", type=int, default=20, help="Number of top TNF matches to return (default: 20)"
)
@click.option(
    "--ani-threshold", type=float, default=80.0, help="Minimum ANI to report (default: 80%)"
)
@click.option("--threads", type=int, default=8, help="Number of threads")
@click.option("--skip-skani", is_flag=True, help="Skip skani ANI calculation (TNF only)")
@click.option("--skip-tnf", is_flag=True, help="Skip TNF distance calculation (ANI only)")
def analyze(
    query, ref, outdir, distance_metric, top_k, ani_threshold, threads, skip_skani, skip_tnf
):
    """Comprehensive genome placement analysis using TNF distance and skani ANI."""
    try:
        if skip_skani and skip_tnf:
            click.echo("\n❌ Error: Cannot skip both skani and TNF analysis", err=True)
            sys.exit(1)

        analyze_genome.analyze_genome(
            query_fna=query,
            ref_fna=ref,
            output_dir=outdir,
            distance_metric=distance_metric,
            top_k=top_k,
            ani_threshold=ani_threshold,
            threads=threads,
            skip_skani=skip_skani,
            skip_tnf=skip_tnf,
        )

    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        import traceback

        traceback.print_exc()
        sys.exit(1)


def main():
    """Entry point."""
    cli()


if __name__ == "__main__":
    main()
