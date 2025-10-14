"""
Canonical 136-dimensional Tetranucleotide Frequency (TNF) calculation.

This module implements the standard reverse-complement collapsed TNF
used in metagenome binning (PMC11289683). All 256 4-mers collapse to
136 canonical representatives (16 self-RC + 240/2 pairs).
"""

import numpy as np
from numba import njit
from numpy.typing import NDArray

# Nucleotide encoding: A=0, C=1, G=2, T=3
_ALPHA = np.frombuffer(b"ACGT", dtype=np.uint8)
_MAP = np.full(256, -1, dtype=np.int16)
_MAP[ord("A")] = 0
_MAP[ord("C")] = 1
_MAP[ord("G")] = 2
_MAP[ord("T")] = 3


def _revcomp_code(x: int) -> int:
    """Reverse complement a 4-mer encoded as 8-bit integer (2 bits per base)."""
    b0 = 3 - (x & 3)
    b1 = 3 - ((x >> 2) & 3)
    b2 = 3 - ((x >> 4) & 3)
    b3 = 3 - ((x >> 6) & 3)
    return (b0 << 6) | (b1 << 4) | (b2 << 2) | b3


# Precompute canonical index for all 256 4-mers (build once at import)
CANON_IDX: NDArray[np.int16] = np.empty(256, dtype=np.int16)
CANON_LIST: list[int] = []  # list of unique canonical 4-mers
seen: dict[int, int] = {}

for x in range(256):
    rc = _revcomp_code(x)
    can = min(x, rc)
    if can not in seen:
        seen[can] = len(CANON_LIST)
        CANON_LIST.append(can)
    CANON_IDX[x] = seen[can]

TNF_DIM = len(CANON_LIST)  # Should be 136


TNFVector = NDArray[np.float32]
ByteArray = NDArray[np.uint8]


@njit(cache=True)
def tnf_136(seq_bytes: ByteArray) -> TNFVector:
    """
    Compute 136-dimensional canonical TNF for a nucleotide sequence.

    Args:
        seq_bytes: Sequence as numpy array of uint8 (ASCII bytes)

    Returns:
        Array of 136 frequencies (sum = 1.0), or zeros if no valid 4-mers

    Notes:
        - Windows containing non-ACGT bases are skipped
        - Reverse-complement collapsed: ACGT and TGCA map to same index
        - Normalized to frequencies (counts / total_valid_windows)
    """
    out: TNFVector = np.zeros(TNF_DIM, dtype=np.float32)
    n = 0
    code = 0
    valid = 0

    for i in range(seq_bytes.size):
        v = _MAP[seq_bytes[i]]
        if v == -1:  # Reset window on ambiguous base
            valid = 0
            continue
        code = ((code << 2) | v) & 0xFF
        if valid < 3:
            valid += 1
            continue
        out[CANON_IDX[code]] += 1.0
        n += 1

    if n > 0:
        out /= n
    return out


def calculate_tnf(sequence: str) -> TNFVector:
    """
    Calculate canonical 136-dim TNF for a DNA sequence string.

    Args:
        sequence: DNA sequence (A, C, G, T, N, etc.)

    Returns:
        NumPy array of 136 frequencies
    """
    seq_upper = sequence.upper().replace("U", "T")
    seq_bytes = np.frombuffer(seq_upper.encode("ascii"), dtype=np.uint8)
    return tnf_136(seq_bytes)


def calculate_tnf_from_contigs(contig_sequences: list[str]) -> TNFVector:
    """
    Calculate genome-level TNF by concatenating all contig sequences.

    Args:
        contig_sequences: List of contig DNA sequences

    Returns:
        NumPy array of 136 frequencies
    """
    concatenated = "".join(contig_sequences)
    return calculate_tnf(concatenated)


def tnf_to_list(tnf_array: TNFVector) -> list[float]:
    """Convert TNF array to list for Parquet storage."""
    return tnf_array.tolist()


def list_to_tnf(tnf_list: list[float]) -> TNFVector:
    """Convert list to TNF array."""
    if len(tnf_list) != TNF_DIM:
        raise ValueError(f"Expected list of length {TNF_DIM}, got {len(tnf_list)}")
    return np.array(tnf_list, dtype=np.float32)
