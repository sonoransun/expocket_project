"""Apply and remove chemical modifications on variants."""

from __future__ import annotations

import numpy as np

from viewer.data.schema import EnrichedVariantDataset, ModificationState
from viewer.encoding.modification_db import MODIFICATIONS_DB, applicable_modifications
from viewer.encoding.property_calculator import compute_modified_properties


class ModificationEngine:
    """Manages chemical modification state for an enriched dataset."""

    def __init__(self, dataset: EnrichedVariantDataset) -> None:
        self._dataset = dataset

    def apply_modification(
        self, variant_id: str, position: int, mod_code: str
    ) -> ModificationState:
        """Apply a modification at a given position."""
        variant = self._dataset.get_variant(variant_id)
        if variant is None:
            raise ValueError(f"Unknown variant: {variant_id}")
        if position < 0 or position >= len(variant.pre_mirna_sequence):
            raise ValueError(f"Position {position} out of range")
        mod = MODIFICATIONS_DB.get(mod_code)
        if mod is None:
            raise ValueError(f"Unknown modification: {mod_code}")
        nt = variant.pre_mirna_sequence[position].upper()
        if nt not in mod.applies_to:
            raise ValueError(f"{mod_code} cannot be applied to {nt}")

        state = self._dataset.get_modification_state(variant_id)
        state.apply(position, mod_code)
        return state

    def remove_modification(
        self, variant_id: str, position: int
    ) -> ModificationState:
        """Remove a modification at a given position."""
        state = self._dataset.get_modification_state(variant_id)
        state.remove(position)
        return state

    def clear_modifications(self, variant_id: str) -> ModificationState:
        """Remove all modifications from a variant."""
        state = self._dataset.get_modification_state(variant_id)
        state.clear()
        return state

    def get_modified_properties(self, variant_id: str) -> np.ndarray:
        """Return (seq_len, 12) property matrix with modifications applied."""
        variant = self._dataset.get_variant(variant_id)
        if variant is None:
            raise ValueError(f"Unknown variant: {variant_id}")
        state = self._dataset.get_modification_state(variant_id)
        return compute_modified_properties(
            variant.pre_mirna_sequence, state.modifications
        )

    def applicable_at_position(
        self, variant_id: str, position: int
    ) -> list[str]:
        """Return modification codes applicable at a given position."""
        variant = self._dataset.get_variant(variant_id)
        if variant is None:
            return []
        if position < 0 or position >= len(variant.pre_mirna_sequence):
            return []
        nt = variant.pre_mirna_sequence[position].upper()
        return [m.code for m in applicable_modifications(nt)]
