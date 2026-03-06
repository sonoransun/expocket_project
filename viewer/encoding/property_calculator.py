"""Compute physicochemical property matrices for datasets."""

from __future__ import annotations

import numpy as np

from viewer.data.schema import VariantDataset
from viewer.encoding.modification_db import MODIFICATIONS_DB
from viewer.encoding.nucleotide_properties import (
    PROPERTY_NAMES,
    encode_sequence,
    get_property,
)


def compute_variant_properties(sequence: str) -> np.ndarray:
    """Encode a single variant's sequence as an (N, 12) property matrix."""
    return encode_sequence(sequence)


def compute_modified_properties(
    sequence: str,
    modifications: dict[int, str],
) -> np.ndarray:
    """Encode a sequence with chemical modifications applied.

    Modifications are applied as property deltas on top of the canonical
    nucleotide properties at each modified position.
    """
    matrix = encode_sequence(sequence)
    for pos, mod_code in modifications.items():
        if pos < 0 or pos >= len(sequence):
            continue
        mod = MODIFICATIONS_DB.get(mod_code)
        if mod is None:
            continue
        # Apply deltas to the canonical property vector at this position
        for prop_name, delta in mod.property_deltas.items():
            if prop_name in PROPERTY_NAMES:
                idx = PROPERTY_NAMES.index(prop_name)
                matrix[pos, idx] += delta
    return matrix


def compute_dataset_property_matrix(dataset: VariantDataset) -> np.ndarray:
    """Compute property matrices for all variants in a dataset.

    Returns an (N_variants, max_seq_len, 12) array.
    """
    if not dataset.variants:
        return np.zeros((0, 0, 12), dtype=np.float64)

    max_len = max(len(v.pre_mirna_sequence) for v in dataset.variants)
    n_variants = len(dataset.variants)
    result = np.zeros((n_variants, max_len, 12), dtype=np.float64)

    for i, variant in enumerate(dataset.variants):
        seq = variant.pre_mirna_sequence
        props = encode_sequence(seq)
        result[i, : len(seq), :] = props

    return result


def compute_summary_features(dataset: VariantDataset) -> np.ndarray:
    """Compute a (N_variants, F) summary feature matrix for SAR analysis.

    Extracts biologically relevant features:
    - Mean properties over full sequence (12 features)
    - Properties at 5' terminal position (12 features)
    - Mean properties at 3' randomized positions (12 features)
    - Properties at cleavage-adjacent positions 19-22 (12 features)
    Total: 48 features per variant.
    """
    if not dataset.variants:
        return np.zeros((0, 48), dtype=np.float64)

    n_variants = len(dataset.variants)
    result = np.zeros((n_variants, 48), dtype=np.float64)

    for i, variant in enumerate(dataset.variants):
        seq = variant.pre_mirna_sequence
        props = encode_sequence(seq)  # (N, 12)

        # Mean over full sequence
        result[i, 0:12] = props.mean(axis=0)

        # 5' terminal position
        if len(seq) > 0:
            result[i, 12:24] = props[0]

        # Mean of last 3 positions (randomized 3' end)
        if len(seq) >= 3:
            result[i, 24:36] = props[-3:].mean(axis=0)

        # Mean of cleavage-adjacent positions (19-22, 0-indexed)
        clv_start = min(19, len(seq) - 1)
        clv_end = min(23, len(seq))
        if clv_end > clv_start:
            result[i, 36:48] = props[clv_start:clv_end].mean(axis=0)

    return result


SUMMARY_FEATURE_NAMES = (
    [f"mean_{p}" for p in PROPERTY_NAMES]
    + [f"5p_{p}" for p in PROPERTY_NAMES]
    + [f"3p_rand_{p}" for p in PROPERTY_NAMES]
    + [f"clv_zone_{p}" for p in PROPERTY_NAMES]
)
