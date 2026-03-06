"""Combinatorial double-replacement screening with pruning."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
from viewer.chemistry.modification_engine import ModificationEngine
from viewer.data.schema import EnrichedVariantDataset
from viewer.encoding.synthesis_db import get_monomer


@dataclass
class DoubleScreenResult:
    """Result for a two-position replacement combination."""

    position_1: int
    position_2: int
    mod_1: str
    mod_2: str
    delta_dc21: float
    delta_dc22: float
    delta_ratio: float
    synergy: float  # actual_double - (single_1_ratio + single_2_ratio)
    synthesis_yield: float


class DoubleReplacementScreener:
    """Screen pairs of modifications across positions for a variant."""

    def __init__(
        self,
        dataset: EnrichedVariantDataset,
        predictor: CleavageSitePredictor,
        engine: ModificationEngine,
    ) -> None:
        self._dataset = dataset
        self._predictor = predictor
        self._engine = engine

    def screen_double(
        self,
        variant_id: str,
        positions: list[int] | None = None,
        mod_codes: list[str] | None = None,
        top_n: int = 50,
    ) -> list[DoubleScreenResult]:
        """Screen all pairs of (position, mod) combinations.

        Uses pruning to keep evaluation count manageable:
        - Default positions: cleavage zone (17-25) + 3' randomized (last 3)
        - Default mods: top 4 by average single-position |delta_ratio|
        """
        variant = self._dataset.get_variant(variant_id)
        if variant is None:
            return []

        seq = variant.pre_mirna_sequence
        seq_len = len(seq)

        if positions is None:
            positions = self._get_relevant_positions(seq_len)

        # Pre-compute singles and cache
        single_cache: dict[tuple[int, str], float] = {}

        if mod_codes is None:
            mod_codes = self._rank_modifications(variant_id, positions)

        # Compute all single-position shifts
        for pos in positions:
            applicable = self._engine.applicable_at_position(variant_id, pos)
            for mod in mod_codes:
                if mod not in applicable:
                    continue
                shifts = self._predictor.predict_shift(variant_id, {pos: mod})
                ratio = shifts.get(21, 0.0) - shifts.get(22, 0.0)
                single_cache[(pos, mod)] = ratio

        # Screen all pairs
        results: list[DoubleScreenResult] = []
        single_keys = list(single_cache.keys())

        for (pos1, mod1), (pos2, mod2) in combinations(single_keys, 2):
            if pos1 == pos2:
                continue

            # Double prediction
            double_shifts = self._predictor.predict_shift(
                variant_id, {pos1: mod1, pos2: mod2}
            )
            d21 = double_shifts.get(21, 0.0)
            d22 = double_shifts.get(22, 0.0)
            double_ratio = d21 - d22

            # Synergy
            single_sum = single_cache[(pos1, mod1)] + single_cache[(pos2, mod2)]
            synergy = double_ratio - single_sum

            # Estimate synthesis yield with both mods
            yield_est = self._estimate_yield(seq, {pos1: mod1, pos2: mod2})

            results.append(DoubleScreenResult(
                position_1=pos1, position_2=pos2,
                mod_1=mod1, mod_2=mod2,
                delta_dc21=d21, delta_dc22=d22,
                delta_ratio=double_ratio,
                synergy=synergy,
                synthesis_yield=yield_est,
            ))

        results.sort(key=lambda r: abs(r.delta_ratio), reverse=True)
        return results[:top_n]

    def compare_single_vs_double(
        self,
        variant_id: str,
        pos1: int,
        mod1: str,
        pos2: int,
        mod2: str,
    ) -> dict:
        """Return detailed comparison of original, single, and double replacements.

        Returns dict with keys: original, single_1, single_2, double, synergy.
        Each maps site -> accuracy delta, except synergy which is a float.
        """
        original = {s: 0.0 for s in [20, 21, 22, 23]}
        single_1 = self._predictor.predict_shift(variant_id, {pos1: mod1})
        single_2 = self._predictor.predict_shift(variant_id, {pos2: mod2})
        double = self._predictor.predict_shift(variant_id, {pos1: mod1, pos2: mod2})

        s1_ratio = single_1.get(21, 0.0) - single_1.get(22, 0.0)
        s2_ratio = single_2.get(21, 0.0) - single_2.get(22, 0.0)
        d_ratio = double.get(21, 0.0) - double.get(22, 0.0)

        return {
            "original": original,
            "single_1": single_1,
            "single_2": single_2,
            "double": double,
            "synergy": d_ratio - (s1_ratio + s2_ratio),
        }

    def _get_relevant_positions(self, seq_len: int) -> list[int]:
        """Cleavage zone (17-25) + 3' randomized (last 3)."""
        positions = list(range(17, min(26, seq_len)))
        for i in range(max(0, seq_len - 3), seq_len):
            if i not in positions:
                positions.append(i)
        return positions

    def _rank_modifications(
        self, variant_id: str, positions: list[int], top_k: int = 4
    ) -> list[str]:
        """Return top mod codes by average |delta_ratio| across positions."""
        from viewer.encoding.modification_db import MODIFICATIONS_DB

        mod_scores: dict[str, float] = {}
        mod_counts: dict[str, int] = {}

        for mod_code in MODIFICATIONS_DB:
            total = 0.0
            count = 0
            for pos in positions[:6]:  # Sample first 6 positions for speed
                applicable = self._engine.applicable_at_position(variant_id, pos)
                if mod_code not in applicable:
                    continue
                ratio = self._predictor.predict_dc_ratio_shift(
                    variant_id, {pos: mod_code}
                )
                total += abs(ratio)
                count += 1
            if count > 0:
                mod_scores[mod_code] = total / count
                mod_counts[mod_code] = count

        ranked = sorted(mod_scores, key=lambda m: mod_scores[m], reverse=True)
        return ranked[:top_k]

    @staticmethod
    def _estimate_yield(
        sequence: str, modifications: dict[int, str]
    ) -> float:
        """Estimate synthesis yield with the given modifications."""
        yield_val = 1.0
        for i in range(len(sequence)):
            nt = sequence[i].upper()
            mod = modifications.get(i)
            monomer = get_monomer(nt, mod)
            yield_val *= monomer.coupling_efficiency
        return yield_val
