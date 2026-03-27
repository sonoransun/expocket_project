"""Predict DICER cleavage shift for edited pre-miRNA sequences."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
    from viewer.data.schema import EnrichedVariantDataset


@dataclass
class EditImpact:
    """DICER cleavage impact of a sequence edit."""
    original_accuracy: dict[int, float]   # {20: ..., 21: ..., 22: ..., 23: ...}
    edited_accuracy: dict[int, float]
    delta: dict[int, float]               # edited - original (positive = more of that site)
    dc21_dc22_ratio_delta: float          # delta[21] - delta[22]
    structural_note: str                  # e.g. "edit in cleavage zone"


class EditImpactPredictor:
    """Predict DICER cleavage accuracy shifts for arbitrary edited sequences.

    Reuses the ridge regression weights already fit by CleavageSitePredictor,
    but accepts any RNA sequence rather than a stored variant ID.
    """

    def __init__(
        self,
        predictor: "CleavageSitePredictor",
        dataset: "EnrichedVariantDataset",
    ) -> None:
        self._predictor = predictor
        self._dataset = dataset

    def predict(
        self,
        variant_id: str,
        modified_rna: str,
    ) -> EditImpact:
        """Compare DICER cleavage prediction for original vs. modified_rna.

        Args:
            variant_id:   Identifies the original (unedited) sequence.
            modified_rna: RNA sequence after the edit.

        Returns EditImpact with per-site deltas.
        """
        variant = self._dataset.get_variant(variant_id)
        pred = self._predictor

        if variant is None or not pred._weights:
            zero = {s: 0.0 for s in [20, 21, 22, 23]}
            return EditImpact(zero, zero, zero, 0.0, "unavailable")

        # Original summary features
        from viewer.data.schema import VariantDataset
        from viewer.encoding.property_calculator import compute_summary_features as _csf
        orig_feat = _csf(VariantDataset(variants=[variant]))[0]  # (48,)

        # Modified summary features
        mod_feat = self._compute_summary_features(modified_rna)

        orig_accuracy: dict[int, float] = {}
        edited_accuracy: dict[int, float] = {}
        delta: dict[int, float] = {}

        for site in [20, 21, 22, 23]:
            w = pred._weights.get(site)
            intercept = pred._intercepts.get(site, 0.0)
            if w is None:
                orig_accuracy[site] = 0.0
                edited_accuracy[site] = 0.0
                delta[site] = 0.0
                continue
            orig_pred = float(((orig_feat - pred._mean) / pred._std) @ w) + intercept
            mod_pred = float(((mod_feat - pred._mean) / pred._std) @ w) + intercept
            orig_accuracy[site] = orig_pred
            edited_accuracy[site] = mod_pred
            delta[site] = mod_pred - orig_pred

        ratio_delta = delta.get(21, 0.0) - delta.get(22, 0.0)

        # Structural note: where is the edit relative to cleavage zone?
        orig_rna = variant.pre_mirna_sequence
        note = _structural_note(orig_rna, modified_rna)

        return EditImpact(
            original_accuracy=orig_accuracy,
            edited_accuracy=edited_accuracy,
            delta=delta,
            dc21_dc22_ratio_delta=ratio_delta,
            structural_note=note,
        )

    @staticmethod
    def _compute_summary_features(rna_seq: str) -> np.ndarray:
        """Replicate the 48-feature summary for an arbitrary RNA sequence."""
        from viewer.encoding.nucleotide_properties import encode_sequence
        props = encode_sequence(rna_seq)  # (N, 12)
        n = props.shape[0]
        result = np.zeros(48, dtype=np.float64)
        if n == 0:
            return result
        result[0:12] = props.mean(axis=0)
        result[12:24] = props[0]
        if n >= 3:
            result[24:36] = props[-3:].mean(axis=0)
        clv_start = min(19, n - 1)
        clv_end = min(23, n)
        if clv_end > clv_start:
            result[36:48] = props[clv_start:clv_end].mean(axis=0)
        return result


def _structural_note(original: str, modified: str) -> str:
    """Describe where the edit falls relative to known DICER features."""
    if len(original) == 0:
        return ""
    # Find first difference
    edit_pos = next(
        (i for i, (a, b) in enumerate(zip(original, modified)) if a != b),
        None,
    )
    if edit_pos is None and len(modified) != len(original):
        edit_pos = min(len(original), len(modified))
    if edit_pos is None:
        return "no change detected"

    if edit_pos < 3:
        return f"edit at pos {edit_pos} — 5' terminus (helicase contact)"
    elif 19 <= edit_pos <= 23:
        return f"edit at pos {edit_pos} — cleavage zone (RNaseIII contact)"
    elif edit_pos >= len(original) - 5:
        return f"edit at pos {edit_pos} — 3' overhang (PAZ contact)"
    elif edit_pos >= len(original) - 3:
        return f"edit at pos {edit_pos} — 3' randomized region"
    else:
        return f"edit at pos {edit_pos} — dsRNA stem"
