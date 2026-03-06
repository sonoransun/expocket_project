"""Data classes for variant, cleavage, and dataset representations."""

from __future__ import annotations

from dataclasses import dataclass, field


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
