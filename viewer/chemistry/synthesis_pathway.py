"""Phosphoramidite synthesis planning for modified oligonucleotides.

Models solid-phase RNA synthesis in the 3'→5' direction with
step-by-step yield tracking and compatibility checking.
"""

from __future__ import annotations

from viewer.data.schema import EnrichedVariantDataset, SynthesisPlan, SynthesisStep
from viewer.encoding.synthesis_db import check_adjacent_compatibility, get_monomer


class SynthesisPlanner:
    """Plans step-by-step oligonucleotide synthesis."""

    def __init__(self, dataset: EnrichedVariantDataset) -> None:
        self._dataset = dataset

    def plan_synthesis(
        self,
        variant_id: str,
        modifications: dict[int, str] | None = None,
        scale_nmol: int = 200,
    ) -> SynthesisPlan:
        """Generate a step-by-step synthesis plan in the 3'→5' direction.

        If modifications is None, uses the current modification state from dataset.
        """
        variant = self._dataset.get_variant(variant_id)
        if variant is None:
            return SynthesisPlan(variant_id=variant_id, sequence="")

        seq = variant.pre_mirna_sequence
        if modifications is None:
            state = self._dataset.get_modification_state(variant_id)
            modifications = dict(state.modifications)

        steps: list[SynthesisStep] = []
        cumulative = 1.0
        warnings: list[str] = []

        # Synthesis proceeds 3' to 5' (reverse of sequence indexing)
        for step_num, seq_pos in enumerate(reversed(range(len(seq)))):
            nt = seq[seq_pos].upper()
            mod = modifications.get(seq_pos)
            monomer = get_monomer(nt, mod)

            cumulative *= monomer.coupling_efficiency

            notes_parts: list[str] = []
            if monomer.requires_extended_coupling:
                notes_parts.append("Extended coupling")

            step = SynthesisStep(
                position=seq_pos,
                nucleotide=nt,
                modification=mod,
                monomer_name=monomer.code,
                coupling_efficiency=monomer.coupling_efficiency,
                cumulative_yield=cumulative,
                deprotection=monomer.deprotection_group,
                cost_factor=monomer.cost_factor,
                notes="; ".join(notes_parts),
            )
            steps.append(step)

        # Check adjacent incompatibilities
        for i in range(len(seq) - 1):
            mod_here = modifications.get(i)
            mod_next = modifications.get(i + 1)
            if mod_here and mod_next:
                if not check_adjacent_compatibility(mod_here, mod_next):
                    warnings.append(
                        f"{mod_here} at pos {i} adjacent to {mod_next} at pos {i+1} "
                        f"— may require modified coupling protocol"
                    )

        total_cost = sum(s.cost_factor for s in steps)

        return SynthesisPlan(
            variant_id=variant_id,
            sequence=seq,
            modifications=modifications,
            steps=steps,
            total_yield=cumulative,
            total_cost_factor=total_cost,
            incompatibilities=warnings,
            scale_nmol=scale_nmol,
        )

    @staticmethod
    def estimate_yield(steps: list[SynthesisStep]) -> float:
        """Product of all coupling efficiencies."""
        result = 1.0
        for s in steps:
            result *= s.coupling_efficiency
        return result

    @staticmethod
    def estimate_cost(steps: list[SynthesisStep]) -> float:
        """Sum of step cost factors."""
        return sum(s.cost_factor for s in steps)
