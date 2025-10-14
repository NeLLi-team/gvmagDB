#!/usr/bin/env python
"""Check if Parquet files have gvclass and annotation columns."""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))


def main() -> None:
    from gvmagdb.catalog import connect

    conn = connect(parquet_glob=str(project_root / "artifacts/parquet/sequences/**/*.parquet"))
    cols = conn.execute("DESCRIBE sequences").df()

    print("Current Parquet Schema:")
    print(f"Total columns: {len(cols)}")
    print()

    # Expected total: 72 (base + gvclass) + 23 (annotations) = 95
    expected_columns = 95

    # Check for gvclass columns
    gvclass_pattern = "gvclass|taxonomy_majority|avgdist|phage_unique|cellular_unique|contigs"
    gvclass_cols = cols[cols["column_name"].str.contains(gvclass_pattern, case=False, na=False)]

    print(f"GVClass columns found: {len(gvclass_cols)}")
    if len(gvclass_cols) > 0:
        print("✅ GVClass columns present:")
        for col in gvclass_cols["column_name"].values:
            print(f"  - {col}")
    else:
        print("❌ NO GVCLASS COLUMNS FOUND")

    print()

    # Check for annotation columns
    annotation_cols = [
        "is_skani_representative",
        "orthogroup_id",
        "is_singleton",
        "emapper_seed_ortholog",
        "emapper_evalue",
        "emapper_score",
        "emapper_eggNOG_OGs",
        "emapper_COG_category",
        "emapper_Description",
        "emapper_Preferred_name",
        "emapper_GOs",
        "emapper_EC",
        "emapper_KEGG_ko",
        "emapper_KEGG_Pathway",
        "emapper_PFAMs",
    ]

    annotation_cols_present = [col for col in annotation_cols if col in cols["column_name"].values]

    print(f"Annotation columns found: {len(annotation_cols_present)}/{len(annotation_cols)}")
    if len(annotation_cols_present) == len(annotation_cols):
        print("✅ All annotation columns present:")
        for col in annotation_cols_present[:5]:  # Show first 5
            print(f"  - {col}")
        print(f"  ... and {len(annotation_cols_present) - 5} more")
    elif len(annotation_cols_present) > 0:
        print("⚠️  PARTIAL annotation columns present:")
        for col in annotation_cols_present:
            print(f"  - {col}")
        missing = set(annotation_cols) - set(annotation_cols_present)
        print("  Missing:")
        for col in missing:
            print(f"  - {col}")
    else:
        print("❌ NO ANNOTATION COLUMNS FOUND")

    print()

    # Final verdict
    if (
        len(cols) == expected_columns
        and len(gvclass_cols) > 0
        and len(annotation_cols_present) == len(annotation_cols)
    ):
        print(f"✅ SCHEMA COMPLETE: {len(cols)} columns (expected {expected_columns})")
        print("✅ RE-INGESTION WITH ANNOTATIONS COMPLETED")
    elif len(cols) == 72 and len(gvclass_cols) > 0:
        print(f"⚠️  SCHEMA INCOMPLETE: {len(cols)} columns (expected {expected_columns})")
        print("⚠️  GVClass present but annotations missing - RE-INGESTION NEEDED")
    elif len(cols) < 72:
        print(f"❌ SCHEMA INCOMPLETE: {len(cols)} columns (expected {expected_columns})")
        print("❌ RE-INGESTION NEEDED")
    else:
        print(f"⚠️  UNEXPECTED SCHEMA: {len(cols)} columns (expected {expected_columns})")

    conn.close()


if __name__ == "__main__":
    main()
