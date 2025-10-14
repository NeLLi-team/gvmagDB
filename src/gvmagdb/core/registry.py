"""
Product registry for tracking derived products.

Tracks DIAMOND databases, skani sketches, HMMER databases, etc.
with versions and checksums for reproducibility.
"""

import hashlib
from datetime import datetime
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pyarrow as pa

# Registry schema
REGISTRY_SCHEMA = pa.schema(
    [
        ("dataset_id", pa.string()),
        ("product", pa.string()),  # 'dmnd', 'skani_db', 'hmm_faa', etc.
        ("uri", pa.string()),  # Path or S3 URI
        ("sha256", pa.string()),  # Checksum
        ("version", pa.string()),  # Tool version
        ("created_at", pa.timestamp("us")),
    ]
)


def compute_sha256(file_path: Path) -> str:
    """Compute SHA256 checksum of a file."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def upsert_product(
    dataset_id: str,
    product: str,
    uri: str,
    version: str,
    registry_path: Path = Path("data/registry/registry.parquet"),
    compute_checksum: bool = True,
) -> None:
    """
    Add or update a product entry in the registry.

    Args:
        dataset_id: Dataset ID
        product: Product type ('dmnd', 'skani_db', etc.)
        uri: Path or S3 URI to product
        version: Tool version (e.g., 'diamond-2.1.9')
        registry_path: Path to registry Parquet file
        compute_checksum: Whether to compute SHA256 checksum
    """
    # Compute checksum if local file
    sha256 = None
    if compute_checksum and Path(uri).exists():
        print(f"  Computing checksum for {uri}...")
        sha256 = compute_sha256(Path(uri))

    # Load existing registry or create new
    if registry_path.exists():
        df = pd.read_parquet(registry_path)
    else:
        df = pd.DataFrame(
            columns=["dataset_id", "product", "uri", "sha256", "version", "created_at"]
        )

    # Remove existing entry for this dataset+product
    df = df[~((df["dataset_id"] == dataset_id) & (df["product"] == product))]

    # Add new entry
    new_entry = pd.DataFrame(
        [
            {
                "dataset_id": dataset_id,
                "product": product,
                "uri": uri,
                "sha256": sha256,
                "version": version,
                "created_at": datetime.now(),
            }
        ]
    )

    df = pd.concat([df, new_entry], ignore_index=True)

    # Write back to Parquet
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(registry_path, schema=REGISTRY_SCHEMA, compression="zstd")

    print(f"✓ Registered: {product} for {dataset_id}")


def get_product(
    dataset_id: str,
    product: str,
    registry_path: Path = Path("data/registry/registry.parquet"),
) -> dict | None:
    """
    Retrieve product entry from registry.

    Args:
        dataset_id: Dataset ID
        product: Product type
        registry_path: Path to registry Parquet

    Returns:
        Dict with product info, or None if not found
    """
    if not registry_path.exists():
        return None

    df = pd.read_parquet(registry_path)
    result = df[(df["dataset_id"] == dataset_id) & (df["product"] == product)]

    if result.empty:
        return None

    entry = result.iloc[0].to_dict()
    return cast(dict[str, Any], entry)


def list_products(
    dataset_id: str | None = None,
    product: str | None = None,
    registry_path: Path = Path("data/registry/registry.parquet"),
) -> pd.DataFrame:
    """
    List products in registry.

    Args:
        dataset_id: Optional filter by dataset_id
        product: Optional filter by product type
        registry_path: Path to registry Parquet

    Returns:
        DataFrame with product entries
    """
    if not registry_path.exists():
        return pd.DataFrame(
            columns=["dataset_id", "product", "uri", "sha256", "version", "created_at"]
        )

    df = pd.read_parquet(registry_path)

    if dataset_id:
        df = df[df["dataset_id"] == dataset_id]

    if product:
        df = df[df["product"] == product]

    return df


def verify_checksum(
    dataset_id: str,
    product: str,
    registry_path: Path = Path("data/registry/registry.parquet"),
) -> bool:
    """
    Verify product checksum matches registry.

    Args:
        dataset_id: Dataset ID
        product: Product type
        registry_path: Path to registry Parquet

    Returns:
        True if checksum matches, False otherwise
    """
    entry = get_product(dataset_id, product, registry_path)
    if not entry:
        print(f"Product not found in registry: {dataset_id}/{product}")
        return False

    uri = entry["uri"]
    expected_sha256 = entry["sha256"]

    if not expected_sha256:
        print(f"No checksum recorded for {dataset_id}/{product}")
        return True  # Can't verify, assume OK

    if not Path(uri).exists():
        print(f"Product file not found: {uri}")
        return False

    actual_sha256 = compute_sha256(Path(uri))

    if actual_sha256 != expected_sha256:
        print(f"❌ Checksum mismatch for {dataset_id}/{product}")
        print(f"   Expected: {expected_sha256}")
        print(f"   Actual:   {actual_sha256}")
        return False

    print(f"✓ Checksum verified for {dataset_id}/{product}")
    return True
