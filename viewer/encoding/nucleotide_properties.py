"""Physicochemical property vectors for canonical RNA nucleotides.

Values derived from:
- Saenger, "Principles of Nucleic Acid Structure" (1984)
- Turner nearest-neighbor thermodynamic parameters
- Bloomfield, Crothers & Tinoco, "Nucleic Acids" (2000)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PROPERTY_NAMES = [
    "molecular_weight",
    "vdw_volume",
    "h_bond_donors",
    "h_bond_acceptors",
    "pka",
    "charge_at_ph7",
    "hydrophobicity_index",
    "sugar_pucker_preference",
    "stacking_energy_5p",
    "stacking_energy_3p",
    "propeller_twist",
    "buckle",
]


@dataclass(frozen=True)
class NucleotidePropertyVector:
    """Physicochemical descriptor vector for a single nucleotide."""

    symbol: str
    display_name: str
    molecular_weight: float       # Da (nucleoside)
    vdw_volume: float             # Angstrom^3
    h_bond_donors: int
    h_bond_acceptors: int
    pka: float                    # most relevant ionizable group
    charge_at_ph7: float          # net charge at physiological pH
    hydrophobicity_index: float   # relative scale (0-1)
    sugar_pucker_preference: float  # 0.0 = C2'-endo, 1.0 = C3'-endo
    stacking_energy_5p: float     # kcal/mol, avg nearest-neighbor 5' context
    stacking_energy_3p: float     # kcal/mol, avg nearest-neighbor 3' context
    propeller_twist: float        # degrees (base-pair geometry)
    buckle: float                 # degrees (base-pair geometry)
    is_purine: bool

    def as_vector(self) -> np.ndarray:
        """Return the 12 numeric properties as a float vector."""
        return np.array(
            [
                self.molecular_weight,
                self.vdw_volume,
                self.h_bond_donors,
                self.h_bond_acceptors,
                self.pka,
                self.charge_at_ph7,
                self.hydrophobicity_index,
                self.sugar_pucker_preference,
                self.stacking_energy_5p,
                self.stacking_energy_3p,
                self.propeller_twist,
                self.buckle,
            ],
            dtype=np.float64,
        )

    @staticmethod
    def vector_size() -> int:
        return 12


# Literature-derived values for canonical RNA nucleotides
CANONICAL_PROPERTIES: dict[str, NucleotidePropertyVector] = {
    "A": NucleotidePropertyVector(
        symbol="A",
        display_name="Adenosine",
        molecular_weight=267.24,
        vdw_volume=135.0,
        h_bond_donors=1,
        h_bond_acceptors=3,
        pka=3.5,
        charge_at_ph7=0.0,
        hydrophobicity_index=0.35,
        sugar_pucker_preference=0.75,  # prefers C3'-endo in RNA
        stacking_energy_5p=-1.70,
        stacking_energy_3p=-1.50,
        propeller_twist=-11.8,
        buckle=0.5,
        is_purine=True,
    ),
    "U": NucleotidePropertyVector(
        symbol="U",
        display_name="Uridine",
        molecular_weight=244.20,
        vdw_volume=112.0,
        h_bond_donors=1,
        h_bond_acceptors=2,
        pka=9.2,
        charge_at_ph7=0.0,
        hydrophobicity_index=0.45,
        sugar_pucker_preference=0.70,
        stacking_energy_5p=-0.90,
        stacking_energy_3p=-1.00,
        propeller_twist=-12.5,
        buckle=-2.4,
        is_purine=False,
    ),
    "G": NucleotidePropertyVector(
        symbol="G",
        display_name="Guanosine",
        molecular_weight=283.24,
        vdw_volume=143.0,
        h_bond_donors=2,
        h_bond_acceptors=4,
        pka=2.1,
        charge_at_ph7=0.0,
        hydrophobicity_index=0.25,
        sugar_pucker_preference=0.80,
        stacking_energy_5p=-2.10,
        stacking_energy_3p=-1.80,
        propeller_twist=-10.4,
        buckle=1.2,
        is_purine=True,
    ),
    "C": NucleotidePropertyVector(
        symbol="C",
        display_name="Cytidine",
        molecular_weight=243.22,
        vdw_volume=108.0,
        h_bond_donors=1,
        h_bond_acceptors=3,
        pka=4.2,
        charge_at_ph7=0.0,
        hydrophobicity_index=0.40,
        sugar_pucker_preference=0.72,
        stacking_energy_5p=-1.30,
        stacking_energy_3p=-1.20,
        propeller_twist=-13.0,
        buckle=-0.8,
        is_purine=False,
    ),
}

# DNA thymine alias (for compatibility with existing data that uses T)
CANONICAL_PROPERTIES["T"] = NucleotidePropertyVector(
    symbol="T",
    display_name="Thymidine (DNA)",
    molecular_weight=242.23,
    vdw_volume=120.0,
    h_bond_donors=1,
    h_bond_acceptors=2,
    pka=9.9,
    charge_at_ph7=0.0,
    hydrophobicity_index=0.50,
    sugar_pucker_preference=0.30,  # prefers C2'-endo in DNA
    stacking_energy_5p=-0.85,
    stacking_energy_3p=-0.95,
    propeller_twist=-12.0,
    buckle=-2.0,
    is_purine=False,
)


def get_property(nt: str) -> NucleotidePropertyVector:
    """Get the property vector for a nucleotide, defaulting to uridine."""
    return CANONICAL_PROPERTIES.get(nt.upper(), CANONICAL_PROPERTIES["U"])


def encode_sequence(sequence: str) -> np.ndarray:
    """Encode a nucleotide sequence as an (N, 12) property matrix."""
    vectors = [get_property(nt).as_vector() for nt in sequence]
    return np.array(vectors, dtype=np.float64)
