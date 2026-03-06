"""Data classes for variant, cleavage, and dataset representations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from viewer.encoding.protein_descriptors import DicerPocketModel


@dataclass
class VariantInfo:
    """Metadata for a single pre-mir-324 variant."""

    variant: str
    group: str  # A, T, G, or C (5' terminal nucleotide)
    randomized_nts: str  # 3-nt randomized 3' end
    pre_mirna_sequence: str  # full ~60 nt sequence
    concrete_struct: str  # dot-bracket secondary structure
    flanking_length_5p: int  # 0 or 1
    new_define_structure_1: str = ""
    new_define_structure_2: str = ""


@dataclass
class CleavageRecord:
    """Cleavage data for one variant at one cleavage site."""

    variant: str
    cleavage_site: int  # 20, 21, 22, or 23
    mean_accuracy: float
    mean_positional_efficiency: float = 0.0
    mean_global_efficiency: float = 0.0
    accuracy_rep1: float = 0.0
    accuracy_rep2: float = 0.0
    accuracy_rep3: float = 0.0


@dataclass
class VariantDataset:
    """Complete dataset ready for visualization."""

    variants: list[VariantInfo] = field(default_factory=list)
    cleavage_data: dict[str, list[CleavageRecord]] = field(default_factory=dict)
    enzyme: str = "hdicer"

    @property
    def variant_ids(self) -> list[str]:
        return [v.variant for v in self.variants]

    def get_variant(self, variant_id: str) -> VariantInfo | None:
        for v in self.variants:
            if v.variant == variant_id:
                return v
        return None

    def get_cleavage(self, variant_id: str, site: int) -> CleavageRecord | None:
        for rec in self.cleavage_data.get(variant_id, []):
            if rec.cleavage_site == site:
                return rec
        return None


@dataclass
class ModificationState:
    """Tracks chemical modifications applied to a single variant."""

    variant: str
    modifications: dict[int, str] = field(default_factory=dict)  # position -> mod code

    def apply(self, position: int, mod_code: str) -> None:
        self.modifications[position] = mod_code

    def remove(self, position: int) -> None:
        self.modifications.pop(position, None)

    def clear(self) -> None:
        self.modifications.clear()


@dataclass
class EnrichedVariantDataset(VariantDataset):
    """VariantDataset extended with physicochemical encodings and modifications."""

    property_matrix: np.ndarray | None = None  # (N_variants, max_seq_len, 12)
    summary_features: np.ndarray | None = None  # (N_variants, 48)
    modification_states: dict[str, ModificationState] = field(default_factory=dict)
    dicer_pocket: DicerPocketModel | None = None

    def get_modification_state(self, variant_id: str) -> ModificationState:
        if variant_id not in self.modification_states:
            self.modification_states[variant_id] = ModificationState(variant=variant_id)
        return self.modification_states[variant_id]


@dataclass
class SynthesisStep:
    """One coupling step in phosphoramidite oligonucleotide synthesis."""

    position: int  # 0-indexed from 3' (synthesis direction)
    nucleotide: str
    modification: str | None
    monomer_name: str
    coupling_efficiency: float
    cumulative_yield: float
    deprotection: str
    cost_factor: float
    notes: str = ""


@dataclass
class SynthesisPlan:
    """Complete synthesis plan for a modified oligonucleotide."""

    variant_id: str
    sequence: str
    modifications: dict[int, str] = field(default_factory=dict)
    steps: list[SynthesisStep] = field(default_factory=list)
    total_yield: float = 0.0
    total_cost_factor: float = 0.0
    incompatibilities: list[str] = field(default_factory=list)
    scale_nmol: int = 200


@dataclass
class ReplacementResult:
    """Result of a single or double nucleotide/modification replacement."""

    positions: list[int] = field(default_factory=list)
    original_bases: list[str] = field(default_factory=list)
    replacements: list[str] = field(default_factory=list)
    shifts: dict[int, float] = field(default_factory=dict)
    dc_ratio_shift: float = 0.0
    synthesis_feasible: bool = True
    estimated_yield: float = 0.0
