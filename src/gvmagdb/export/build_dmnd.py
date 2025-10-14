"""
DIAMOND database builder.

Exports proteins to FASTA and builds DIAMOND .dmnd database.
"""

import subprocess
from pathlib import Path

from gvmagdb.core.registry import upsert_product
from gvmagdb.export.export_faa import export_aa_fasta


def build_diamond_db(
    dataset_id: str | None = None,
    dmnd_output: Path = Path("artifacts/products/dmnd/all_proteins.dmnd"),
    threads: int = 8,
) -> Path:
    """
    Build DIAMOND database from protein sequences.

    Args:
        dataset_id: Optional filter by dataset (None = all proteins)
        dmnd_output: Output .dmnd file path
        threads: Number of threads

    Returns:
        Path to created .dmnd file
    """
    print("\n" + "=" * 70)
    print("Building DIAMOND Database")
    print("=" * 70)

    # Export proteins to temp FASTA
    faa_temp = dmnd_output.parent / "proteins_temp.faa"
    faa_temp.parent.mkdir(parents=True, exist_ok=True)

    print(f"\n📤 Exporting proteins to {faa_temp}...")
    count = export_aa_fasta(faa_temp, dataset_id=dataset_id)

    if count == 0:
        raise ValueError("No proteins found to export!")

    # Build DIAMOND database
    print("\n🔨 Building DIAMOND database...")
    dmnd_base = str(dmnd_output.with_suffix(""))  # diamond adds .dmnd

    cmd = [
        "diamond",
        "makedb",
        "--in",
        str(faa_temp),
        "--db",
        dmnd_base,
        "--threads",
        str(threads),
    ]

    print(f"   Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Clean up temp FASTA
    faa_temp.unlink()

    # Register in registry
    dataset_label = dataset_id if dataset_id else "all"
    upsert_product(
        dataset_id=dataset_label,
        product="dmnd",
        uri=str(dmnd_output),
        version=get_diamond_version(),
    )

    print(f"\n✓ DIAMOND database created: {dmnd_output}")
    print(f"✓ Indexed {count} protein sequences")
    print("=" * 70 + "\n")

    return dmnd_output


def get_diamond_version() -> str:
    """Get DIAMOND version string."""
    try:
        result = subprocess.run(
            ["diamond", "version"],
            capture_output=True,
            text=True,
            check=True,
        )
        # Parse version from output
        for line in result.stdout.split("\n"):
            if "diamond version" in line.lower():
                return line.strip()
        return "diamond-unknown"
    except Exception:
        return "diamond-unknown"
