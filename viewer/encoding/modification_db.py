"""Chemical modification database for non-standard nucleotides.

Each modification is defined as property deltas relative to the parent
canonical nucleotide. Modified properties = canonical + deltas.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ChemicalModification:
    """A non-standard nucleotide chemical modification."""

    code: str
    full_name: str
    applies_to: tuple[str, ...]  # canonical nts this can modify
    property_deltas: dict[str, float]  # changes relative to parent
    stability_effect: float  # Tm shift per modification (degrees)
    dicer_recognition_effect: float  # -1 to +1 (negative = reduced)
    nuclease_resistance: float  # 0-1 scale
    display_shape: str  # 3D geometry hint
    display_color: tuple[float, float, float, float]  # RGBA
    description: str


MODIFICATIONS_DB: dict[str, ChemicalModification] = {
    "2OMe": ChemicalModification(
        code="2OMe",
        full_name="2'-O-Methylation",
        applies_to=("A", "U", "G", "C"),
        property_deltas={
            "molecular_weight": 14.0,
            "vdw_volume": 15.0,
            "sugar_pucker_preference": 0.15,
            "hydrophobicity_index": 0.12,
        },
        stability_effect=1.5,
        dicer_recognition_effect=-0.3,
        nuclease_resistance=0.7,
        display_shape="cone",
        display_color=(0.2, 0.8, 0.9, 0.85),
        description=(
            "Methyl group on 2'-OH of ribose. Increases duplex stability "
            "and nuclease resistance. Reduces DICER recognition when near "
            "the cleavage site."
        ),
    ),
    "LNA": ChemicalModification(
        code="LNA",
        full_name="Locked Nucleic Acid",
        applies_to=("A", "U", "G", "C"),
        property_deltas={
            "molecular_weight": 12.0,
            "vdw_volume": 18.0,
            "sugar_pucker_preference": 0.25,
            "stacking_energy_5p": -0.50,
            "stacking_energy_3p": -0.40,
        },
        stability_effect=4.0,
        dicer_recognition_effect=-0.6,
        nuclease_resistance=0.9,
        display_shape="cube",
        display_color=(0.9, 0.3, 0.1, 0.85),
        description=(
            "Methylene bridge locking C3'-endo sugar pucker. Very high "
            "binding affinity (+4°C Tm per mod). Strongly inhibits DICER "
            "cleavage when positioned near the active site."
        ),
    ),
    "PSI": ChemicalModification(
        code="PSI",
        full_name="Pseudouridine (Ψ)",
        applies_to=("U",),
        property_deltas={
            "h_bond_donors": 1,
            "stacking_energy_5p": -0.30,
            "hydrophobicity_index": -0.05,
        },
        stability_effect=1.0,
        dicer_recognition_effect=-0.1,
        nuclease_resistance=0.3,
        display_shape="dodecahedron",
        display_color=(0.7, 0.5, 0.9, 0.85),
        description=(
            "C-glycosidic isomer of uridine. Gains an extra H-bond donor. "
            "Slightly increases thermal stability with minimal impact on "
            "DICER processing."
        ),
    ),
    "m6A": ChemicalModification(
        code="m6A",
        full_name="N6-Methyladenosine",
        applies_to=("A",),
        property_deltas={
            "molecular_weight": 14.0,
            "h_bond_donors": -1,
            "stacking_energy_5p": 0.30,
            "hydrophobicity_index": 0.10,
        },
        stability_effect=-1.5,
        dicer_recognition_effect=-0.4,
        nuclease_resistance=0.1,
        display_shape="octahedron",
        display_color=(0.9, 0.8, 0.2, 0.85),
        description=(
            "Methyl group on the N6 amino group of adenosine. Weakens "
            "base pairing by disrupting Watson-Crick hydrogen bonding. "
            "Destabilizes duplex and reduces DICER recognition."
        ),
    ),
    "2F": ChemicalModification(
        code="2F",
        full_name="2'-Fluoro",
        applies_to=("A", "U", "G", "C"),
        property_deltas={
            "molecular_weight": 2.0,
            "vdw_volume": -2.0,
            "sugar_pucker_preference": 0.20,
            "hydrophobicity_index": 0.08,
        },
        stability_effect=2.0,
        dicer_recognition_effect=-0.2,
        nuclease_resistance=0.8,
        display_shape="cone",
        display_color=(0.3, 0.9, 0.4, 0.85),
        description=(
            "Fluorine replacing 2'-OH. Locks sugar in C3'-endo pucker. "
            "High thermal stability and nuclease resistance with moderate "
            "DICER compatibility."
        ),
    ),
    "PS": ChemicalModification(
        code="PS",
        full_name="Phosphorothioate",
        applies_to=("A", "U", "G", "C"),
        property_deltas={
            "molecular_weight": 16.0,
            "charge_at_ph7": -0.1,
            "hydrophobicity_index": 0.15,
        },
        stability_effect=-0.5,
        dicer_recognition_effect=-0.15,
        nuclease_resistance=0.85,
        display_shape="cone",
        display_color=(1.0, 0.6, 0.0, 0.85),
        description=(
            "Sulfur replacing non-bridging oxygen in the phosphodiester "
            "backbone. Slightly destabilizes duplex but confers strong "
            "nuclease resistance. Backbone modification, not base."
        ),
    ),
    "s4U": ChemicalModification(
        code="s4U",
        full_name="4-Thiouridine",
        applies_to=("U",),
        property_deltas={
            "molecular_weight": 16.0,
            "vdw_volume": 8.0,
            "stacking_energy_5p": -0.20,
            "hydrophobicity_index": 0.20,
        },
        stability_effect=-1.0,
        dicer_recognition_effect=-0.2,
        nuclease_resistance=0.2,
        display_shape="octahedron",
        display_color=(0.8, 0.7, 0.1, 0.85),
        description=(
            "Sulfur replacing the C4 oxygen in uracil. UV-crosslinkable. "
            "Slightly destabilizes base pairing due to larger atomic radius "
            "of sulfur. Useful as a photo-crosslinking probe."
        ),
    ),
    "m5C": ChemicalModification(
        code="m5C",
        full_name="5-Methylcytosine",
        applies_to=("C",),
        property_deltas={
            "molecular_weight": 14.0,
            "stacking_energy_5p": -0.25,
            "hydrophobicity_index": 0.08,
        },
        stability_effect=0.5,
        dicer_recognition_effect=-0.05,
        nuclease_resistance=0.1,
        display_shape="dodecahedron",
        display_color=(0.5, 0.8, 0.7, 0.85),
        description=(
            "Methyl group at the C5 position of cytosine. Slightly "
            "enhances stacking. Common epigenetic mark in DNA, also found "
            "in RNA. Minimal effect on DICER processing."
        ),
    ),
    "INO": ChemicalModification(
        code="INO",
        full_name="Inosine",
        applies_to=("A",),
        property_deltas={
            "molecular_weight": -1.0,
            "h_bond_donors": -1,
            "h_bond_acceptors": 0,
            "stacking_energy_5p": 0.40,
        },
        stability_effect=-2.0,
        dicer_recognition_effect=-0.35,
        nuclease_resistance=0.05,
        display_shape="octahedron",
        display_color=(0.6, 0.4, 0.8, 0.85),
        description=(
            "Deaminated adenosine (A-to-I editing product). Pairs with "
            "cytidine instead of uridine. Destabilizes A-U pairs and "
            "can alter DICER cleavage site selection."
        ),
    ),
}


def get_modification(code: str) -> ChemicalModification | None:
    """Look up a modification by its code."""
    return MODIFICATIONS_DB.get(code)


def applicable_modifications(nucleotide: str) -> list[ChemicalModification]:
    """Return all modifications applicable to a given canonical nucleotide."""
    nt = nucleotide.upper()
    return [m for m in MODIFICATIONS_DB.values() if nt in m.applies_to]
