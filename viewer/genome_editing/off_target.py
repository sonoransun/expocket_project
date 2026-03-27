"""Off-target scoring against all variants in the dataset."""

from __future__ import annotations

from typing import TYPE_CHECKING

from viewer.genome_editing.sequence_utils import (
    get_dual_strand,
    find_pam_sites,
    count_mismatches_weighted,
)

if TYPE_CHECKING:
    from viewer.data.schema import EnrichedVariantDataset

# Position weights for 20-nt guide: seed (positions 0-7, PAM-proximal) weight 2.0,
# rest 0.5.  Position 0 = PAM-adjacent nt.
_DEFAULT_WEIGHTS: list[float] = [2.0 if i < 8 else 0.5 for i in range(20)]


def _extract_spacers_from_dna(
    sense_dna: str,
    pam_pattern: str | None,
    pam_position: str,
    spacer_len: int,
) -> list[tuple[str, str, int]]:
    """Find all potential target spacers in a DNA sequence.

    Returns list of (spacer_seq, strand, start_pos).
    """
    if pam_pattern is None or spacer_len <= 0:
        return []
    matches = find_pam_sites(sense_dna, pam_pattern, pam_position, spacer_len)
    return [(m.spacer_seq, m.strand, m.spacer_start) for m in matches]


def find_off_targets_in_dataset(
    guide_spacer: str,
    pam_pattern: str | None,
    pam_position: str,
    spacer_len: int,
    dataset: "EnrichedVariantDataset",
    exclude_variant_id: str = "",
    threshold: float = 0.8,
) -> list[dict]:
    """Scan all variants for PAM + spacer with ≤ threshold mismatch score.

    threshold is the maximum weighted mismatch score to call something an off-target.
    Lower threshold = stricter (fewer off-targets flagged).

    Returns list of dicts with keys: variant_id, spacer, strand, start, mismatches, weighted_score.
    """
    if pam_pattern is None or spacer_len <= 0:
        return []

    results: list[dict] = []
    guide = guide_spacer.upper()
    n = len(guide)
    weights = [2.0 if i < 8 else 0.5 for i in range(n)]

    for variant in dataset.variants:
        vid = variant.variant
        if vid == exclude_variant_id:
            continue
        sense_dna, _ = get_dual_strand(variant.pre_mirna_sequence)
        spacers = _extract_spacers_from_dna(
            sense_dna, pam_pattern, pam_position, spacer_len
        )
        for spacer, strand, start in spacers:
            ws = count_mismatches_weighted(guide, spacer.upper(), weights)
            if ws <= threshold:
                mm = sum(1 for a, b in zip(guide, spacer.upper()) if a != b)
                results.append({
                    "variant_id": vid,
                    "spacer": spacer,
                    "strand": strand,
                    "start": start,
                    "mismatches": mm,
                    "weighted_score": ws,
                })

    results.sort(key=lambda x: x["weighted_score"])
    return results


def compute_specificity(off_targets: list[dict]) -> float:
    """Return specificity score 0-1 (1 = no off-targets)."""
    if not off_targets:
        return 1.0
    # Penalise each off-target by its similarity
    penalty = sum(max(0.0, 1.0 - r["weighted_score"] / 8.0) for r in off_targets)
    return max(0.0, 1.0 - penalty / (1.0 + len(off_targets)))


class OffTargetScorer:
    """Stateless façade for off-target analysis."""

    def score(
        self,
        guide_spacer: str,
        pam_pattern: str | None,
        pam_position: str,
        spacer_len: int,
        dataset: "EnrichedVariantDataset",
        exclude_variant_id: str = "",
    ) -> tuple[float, list[dict]]:
        """Return (specificity_score, off_target_list)."""
        ots = find_off_targets_in_dataset(
            guide_spacer, pam_pattern, pam_position, spacer_len,
            dataset, exclude_variant_id,
        )
        return compute_specificity(ots), ots
