"""Screen all position x modification combos and rank by predicted impact."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
from viewer.chemistry.modification_engine import ModificationEngine
from viewer.data.schema import EnrichedVariantDataset


@dataclass
class ScreenResult:
    """Result for a single position x modification combo."""

    position: int
    modification: str
    delta_dc21: float
    delta_dc22: float
    delta_ratio: float  # delta_dc21 - delta_dc22


class VirtualScreener:
    """Screen modifications across all positions for a variant."""

    def __init__(
        self,
        dataset: EnrichedVariantDataset,
        engine: ModificationEngine,
        predictor: CleavageSitePredictor,
    ) -> None:
        self._dataset = dataset
        self._engine = engine
        self._predictor = predictor

    def screen_variant(self, variant_id: str) -> list[ScreenResult]:
        """Screen all applicable modifications at all positions.

        Returns results sorted by |delta_ratio| descending.
        """
        variant = self._dataset.get_variant(variant_id)
        if variant is None:
            return []

        results: list[ScreenResult] = []
        seq = variant.pre_mirna_sequence

        for pos in range(len(seq)):
            applicable = self._engine.applicable_at_position(variant_id, pos)
            for mod_code in applicable:
                shifts = self._predictor.predict_shift(
                    variant_id, {pos: mod_code}
                )
                d21 = shifts.get(21, 0.0)
                d22 = shifts.get(22, 0.0)
                results.append(ScreenResult(
                    position=pos,
                    modification=mod_code,
                    delta_dc21=d21,
                    delta_dc22=d22,
                    delta_ratio=d21 - d22,
                ))

        results.sort(key=lambda r: abs(r.delta_ratio), reverse=True)
        return results

    def rank_by_dc_ratio_shift(
        self,
        variant_id: str,
        direction: str = "increase_dc21",
        top_n: int = 20,
    ) -> list[ScreenResult]:
        """Screen and return top results favoring a direction.

        direction: 'increase_dc21' (positive delta_ratio) or
                   'increase_dc22' (negative delta_ratio).
        """
        all_results = self.screen_variant(variant_id)
        if direction == "increase_dc21":
            all_results.sort(key=lambda r: r.delta_ratio, reverse=True)
        else:
            all_results.sort(key=lambda r: r.delta_ratio)
        return all_results[:top_n]
