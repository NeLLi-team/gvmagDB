#!/usr/bin/env python
"""Test annotation queries on the database."""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))


def main() -> None:
    from gvmagdb.catalog import connect

    conn = connect(parquet_glob=str(project_root / "artifacts/parquet/sequences/**/*.parquet"))

    print("=" * 70)
    print("Annotation Integration Test")
    print("=" * 70)

    # Test 1: Count representative vs non-representative proteins
    print("\n1. Representative genome distribution:")
    query1 = """
SELECT
    is_skani_representative,
    seq_type,
    COUNT(*) as count
FROM sequences
GROUP BY is_skani_representative, seq_type
ORDER BY is_skani_representative DESC, seq_type
"""
    result1 = conn.execute(query1).df()
    print(result1.to_string(index=False))

    # Test 2: Count proteins with orthogroup assignments
    print("\n2. Orthogroup assignment coverage:")
    query2 = """
SELECT
    COUNT(*) as total_proteins,
    COUNT(orthogroup_id) as with_orthogroup,
    SUM(CASE WHEN is_singleton THEN 1 ELSE 0 END) as singletons,
    SUM(CASE WHEN orthogroup_id IS NULL AND is_singleton IS NULL THEN 1 ELSE 0 END) as no_assignment
FROM sequences
WHERE seq_type = 'AA'
"""
    result2 = conn.execute(query2).df()
    print(result2.to_string(index=False))

    # Test 3: Count proteins with eggNOG annotations
    print("\n3. eggNOG annotation coverage:")
    query3 = """
SELECT
    COUNT(*) as total_proteins,
    COUNT(emapper_Description) as with_description,
    COUNT(emapper_COG_category) as with_COG,
    COUNT(emapper_KEGG_ko) as with_KEGG,
    COUNT(emapper_PFAMs) as with_PFAM
FROM sequences
WHERE seq_type = 'AA'
"""
    result3 = conn.execute(query3).df()
    print(result3.to_string(index=False))

    # Test 4: Sample annotated proteins
    print("\n4. Sample annotated proteins (first 5):")
    query4 = """
SELECT
    gid,
    dataset_id,
    is_skani_representative,
    orthogroup_id,
    is_singleton,
    emapper_Description,
    emapper_COG_category
FROM sequences
WHERE seq_type = 'AA'
    AND is_skani_representative = true
    AND emapper_Description IS NOT NULL
LIMIT 5
"""
    result4 = conn.execute(query4).df()
    print(result4.to_string(index=False))

    # Test 5: Check for specific functional categories
    print("\n5. Top 10 COG categories:")
    query5 = """
SELECT
    emapper_COG_category,
    COUNT(*) as count
FROM sequences
WHERE seq_type = 'AA' AND emapper_COG_category IS NOT NULL
GROUP BY emapper_COG_category
ORDER BY count DESC
LIMIT 10
"""
    result5 = conn.execute(query5).df()
    print(result5.to_string(index=False))

    print("\n" + "=" * 70)
    print("✅ Annotation integration test complete!")
    print("=" * 70)

    conn.close()


if __name__ == "__main__":
    main()
