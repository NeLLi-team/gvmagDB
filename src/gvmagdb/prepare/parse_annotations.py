#!/usr/bin/env python
"""
Parse protein annotation files for gvmagDB ingestion.

This module provides parsers for:
- ProteinOrtho orthogroup assignments
- ProteinOrtho singleton proteins
- eggNOG-mapper functional annotations
"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_proteinortho(proteinortho_tsv: Path) -> dict[str, str]:
    """
    Parse ProteinOrtho TSV to extract orthogroup assignments.

    Args:
        proteinortho_tsv: Path to gvmags_db_trial4.proteinortho.tsv

    Returns:
        Dictionary mapping protein_id → orthogroup_id (row number as string)

    Format:
        # Species	Genes	Alg.-Conn.	AbALV.faa	Ac22.faa	...
        3436	3478	0.101	*	*	*	...
        3085	3136	0.109	*	Ac22|Ac22_1_315	*	...
    """
    protein_to_orthogroup = {}

    logger.info(f"Parsing ProteinOrtho file: {proteinortho_tsv}")

    with open(proteinortho_tsv) as f:
        # Read header to identify genome columns (skip first 3 columns: Species, Genes, Alg.-Conn.)
        f.readline()
        genome_col_start = 3

        # Process each orthogroup (row)
        for row_num, line in enumerate(f, start=1):
            fields = line.strip().split("\t")
            orthogroup_id = str(row_num)

            # Check each genome column for protein members
            for col_idx in range(genome_col_start, len(fields)):
                cell = fields[col_idx]
                if cell != "*":
                    # Cell format: genome|protein_id
                    protein_id = cell
                    protein_to_orthogroup[protein_id] = orthogroup_id

    logger.info(f"Parsed {len(protein_to_orthogroup):,} proteins in orthogroups")
    return protein_to_orthogroup


def parse_singletons(singletons_tsv: Path) -> set[str]:
    """
    Parse ProteinOrtho singletons TSV to extract singleton proteins.

    Args:
        singletons_tsv: Path to gvmags_db_trial4.singletons.tsv

    Returns:
        Set of protein_ids that are singletons

    Format:
        # Species	Genes	Alg.-Conn.	AbALV.faa	Ac22.faa	...
        1	1	0	AbALV|LC506465__1_21	*	*	...
        1	1	0	AbALV|LC506465__1_23	*	*	...
    """
    singletons = set()

    logger.info(f"Parsing singletons file: {singletons_tsv}")

    with open(singletons_tsv) as f:
        # Skip header
        next(f, None)

        # Process each singleton (each row has exactly one protein)
        for line_num, line in enumerate(f, start=1):
            fields = line.strip().split("\t")

            # Find the non-'*' cell (should be exactly one per row)
            for cell in fields[3:]:  # Skip first 3 columns
                if cell != "*":
                    protein_id = cell
                    singletons.add(protein_id)
                    break

            # Progress indicator for large file
            if line_num % 1_000_000 == 0:
                logger.info(f"  Processed {line_num:,} singleton rows...")

    logger.info(f"Parsed {len(singletons):,} singleton proteins")
    return singletons


def parse_eggnog(eggnog_annotations: Path) -> dict[str, dict[str, object | None]]:
    """
    Parse eggNOG-mapper annotations TSV.

    Args:
        eggnog_annotations: Path to emapper_out.emapper.annotations

    Returns:
        Dictionary mapping protein_id → annotation_dict

    Annotation dict contains these keys:
        seed_ortholog, evalue, score, eggNOG_OGs, max_annot_lvl,
        COG_category, Description, Preferred_name, GOs, EC,
        KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass,
        BRITE, KEGG_TC, CAZy, BiGG_Reaction, PFAMs

    Format:
        #query	seed_ortholog	evalue	score	...
        GVMAG-M-3300023184-120|Ga0214919_10000160_3	...	1e-50	200	...
    """
    protein_annotations: dict[str, dict[str, object | None]] = {}

    logger.info(f"Parsing eggNOG annotations: {eggnog_annotations}")

    # Column names (excluding 'query' which is the protein_id)
    annotation_columns = [
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

    with open(eggnog_annotations) as f:
        header: list[str] | None = None
        col_indices: dict[str, int] | None = None

        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("#query"):
                    header = line.lstrip("#").split("\t")
                    col_indices = {col: idx for idx, col in enumerate(header)}
                continue

            if header is None:
                header = ["query"] + annotation_columns
                col_indices = {col: idx for idx, col in enumerate(header)}

            assert col_indices is not None  # for type checkers

            fields = line.split("\t")
            protein_id = fields[col_indices["query"]] if "query" in col_indices else fields[0]

            annotations: dict[str, object | None] = {}
            for col in annotation_columns:
                idx = col_indices.get(col)
                if idx is None or idx >= len(fields):
                    annotations[col] = None
                    continue

                value = fields[idx]
                if not value or value == "-":
                    annotations[col] = None
                elif col in {"evalue", "score"}:
                    try:
                        annotations[col] = float(value)
                    except (ValueError, TypeError):
                        annotations[col] = None
                else:
                    annotations[col] = value

            protein_annotations[protein_id] = annotations

            if line_num % 500_000 == 0:
                logger.info("  Processed %s annotation rows...", f"{line_num:,}")

    logger.info(f"Parsed annotations for {len(protein_annotations):,} proteins")
    return protein_annotations


def load_all_annotations(proteinortho_tsv: Path, singletons_tsv: Path, eggnog_tsv: Path) -> tuple[
    dict[str, str],
    set[str],
    dict[str, dict[str, object | None]],
]:
    """
    Load all annotation files at once.

    Returns:
        Tuple of (orthogroup_dict, singleton_set, eggnog_dict)
    """
    logger.info("Loading all protein annotations...")

    orthogroup_dict = parse_proteinortho(proteinortho_tsv)
    singleton_set = parse_singletons(singletons_tsv)
    eggnog_dict = parse_eggnog(eggnog_tsv)

    logger.info("All annotations loaded successfully")
    return orthogroup_dict, singleton_set, eggnog_dict
