"""Amino acid property vectors and DICER pocket model.

Amino acid descriptors from Kyte-Doolittle hydrophobicity, Zamyatnin volumes,
and flexibility indices. DICER pocket contacts based on PDB:5ZAL
(human DICER-pre-miRNA complex).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class AminoAcidProperties:
    """Physicochemical descriptor for a single amino acid residue."""

    one_letter: str
    three_letter: str
    hydrophobicity: float  # Kyte-Doolittle scale (-4.5 to 4.5)
    charge_ph7: float
    volume: float  # van der Waals volume (Angstrom^3)
    flexibility: float  # normalized B-factor (0-1)
    h_bond_capacity: float  # combined donors + acceptors (0-1 scale)

    def as_vector(self) -> np.ndarray:
        return np.array(
            [self.hydrophobicity, self.charge_ph7, self.volume,
             self.flexibility, self.h_bond_capacity],
            dtype=np.float64,
        )


AA_PROPERTIES: dict[str, AminoAcidProperties] = {
    "A": AminoAcidProperties("A", "Ala", 1.8, 0.0, 88.6, 0.36, 0.1),
    "R": AminoAcidProperties("R", "Arg", -4.5, 1.0, 173.4, 0.53, 0.9),
    "N": AminoAcidProperties("N", "Asn", -3.5, 0.0, 114.1, 0.46, 0.7),
    "D": AminoAcidProperties("D", "Asp", -3.5, -1.0, 111.1, 0.51, 0.6),
    "C": AminoAcidProperties("C", "Cys", 2.5, 0.0, 108.5, 0.35, 0.2),
    "E": AminoAcidProperties("E", "Glu", -3.5, -1.0, 138.4, 0.50, 0.6),
    "Q": AminoAcidProperties("Q", "Gln", -3.5, 0.0, 143.8, 0.49, 0.7),
    "G": AminoAcidProperties("G", "Gly", -0.4, 0.0, 60.1, 0.54, 0.1),
    "H": AminoAcidProperties("H", "His", -3.2, 0.1, 153.2, 0.41, 0.6),
    "I": AminoAcidProperties("I", "Ile", 4.5, 0.0, 166.7, 0.28, 0.1),
    "L": AminoAcidProperties("L", "Leu", 3.8, 0.0, 166.7, 0.32, 0.1),
    "K": AminoAcidProperties("K", "Lys", -3.9, 1.0, 168.6, 0.47, 0.7),
    "M": AminoAcidProperties("M", "Met", 1.9, 0.0, 162.9, 0.42, 0.2),
    "F": AminoAcidProperties("F", "Phe", 2.8, 0.0, 189.9, 0.31, 0.1),
    "P": AminoAcidProperties("P", "Pro", -1.6, 0.0, 112.7, 0.51, 0.1),
    "S": AminoAcidProperties("S", "Ser", -0.8, 0.0, 89.0, 0.51, 0.5),
    "T": AminoAcidProperties("T", "Thr", -0.7, 0.0, 116.1, 0.44, 0.5),
    "W": AminoAcidProperties("W", "Trp", -0.9, 0.0, 227.8, 0.27, 0.3),
    "Y": AminoAcidProperties("Y", "Tyr", -1.3, 0.0, 193.6, 0.33, 0.4),
    "V": AminoAcidProperties("V", "Val", 4.2, 0.0, 140.0, 0.30, 0.1),
}


@dataclass
class DicerPocketResidue:
    """A single residue in the DICER RNA-binding pocket."""

    residue_id: int
    amino_acid: str  # one-letter code
    domain: str  # PAZ, RNaseIIIa, RNaseIIIb, helicase, platform
    rna_contact_positions: list[int] = field(default_factory=list)
    interaction_type: str = "vdw"  # vdw, hbond, electrostatic, stacking
    interaction_strength: float = 0.5  # 0-1 scale

    @property
    def properties(self) -> AminoAcidProperties:
        return AA_PROPERTIES.get(self.amino_acid, AA_PROPERTIES["G"])


@dataclass
class DicerPocketModel:
    """Simplified model of the DICER enzyme pocket."""

    residues: list[DicerPocketResidue] = field(default_factory=list)
    enzyme: str = "hdicer"

    def contact_matrix(self, seq_length: int) -> np.ndarray:
        """Return (n_residues, seq_length) interaction strength matrix."""
        mat = np.zeros((len(self.residues), seq_length), dtype=np.float64)
        for i, res in enumerate(self.residues):
            for pos in res.rna_contact_positions:
                if 0 <= pos < seq_length:
                    mat[i, pos] = res.interaction_strength
        return mat

    def domain_residues(self, domain: str) -> list[DicerPocketResidue]:
        return [r for r in self.residues if r.domain == domain]

    @property
    def domains(self) -> list[str]:
        return sorted(set(r.domain for r in self.residues))


def build_mock_dicer_pocket(enzyme: str = "hdicer") -> DicerPocketModel:
    """Build a biologically reasonable mock DICER pocket model.

    Based on structural data from PDB:5ZAL (human DICER-pre-miRNA complex).
    The PAZ domain grips the 3' 2-nt overhang (positions ~58-62 in a 63-nt hairpin).
    RNaseIIIa cleaves the 3p strand (~position 20-22).
    RNaseIIIb cleaves the 5p strand (~position 41-43, counting from 5').
    """
    residues = [
        # PAZ domain — recognizes 3' overhang
        DicerPocketResidue(227, "Y", "PAZ", [60, 61], "stacking", 0.8),
        DicerPocketResidue(230, "R", "PAZ", [60, 61, 62], "electrostatic", 0.9),
        DicerPocketResidue(233, "H", "PAZ", [59, 60], "hbond", 0.7),
        DicerPocketResidue(252, "F", "PAZ", [61, 62], "stacking", 0.8),
        DicerPocketResidue(275, "K", "PAZ", [58, 59], "electrostatic", 0.6),
        DicerPocketResidue(278, "R", "PAZ", [57, 58], "electrostatic", 0.7),
        DicerPocketResidue(295, "Y", "PAZ", [59], "stacking", 0.5),
        DicerPocketResidue(310, "N", "PAZ", [60], "hbond", 0.6),
        # Platform domain — measures dsRNA
        DicerPocketResidue(342, "R", "platform", [5, 6], "electrostatic", 0.5),
        DicerPocketResidue(345, "K", "platform", [4, 5], "electrostatic", 0.5),
        DicerPocketResidue(360, "S", "platform", [6, 7], "hbond", 0.4),
        DicerPocketResidue(378, "R", "platform", [3, 4], "electrostatic", 0.6),
        # RNaseIIIa — cleaves 3p strand (near positions 20-22)
        DicerPocketResidue(1316, "E", "RNaseIIIa", [19, 20], "electrostatic", 0.9),
        DicerPocketResidue(1320, "D", "RNaseIIIa", [20, 21], "electrostatic", 0.95),
        DicerPocketResidue(1324, "D", "RNaseIIIa", [21], "electrostatic", 0.9),
        DicerPocketResidue(1328, "E", "RNaseIIIa", [20, 21, 22], "electrostatic", 0.85),
        DicerPocketResidue(1335, "R", "RNaseIIIa", [18, 19], "electrostatic", 0.6),
        DicerPocketResidue(1340, "K", "RNaseIIIa", [22, 23], "electrostatic", 0.5),
        DicerPocketResidue(1345, "N", "RNaseIIIa", [19, 20], "hbond", 0.7),
        # RNaseIIIb — cleaves 5p strand (near positions 40-43)
        DicerPocketResidue(1636, "E", "RNaseIIIb", [40, 41], "electrostatic", 0.9),
        DicerPocketResidue(1640, "D", "RNaseIIIb", [41, 42], "electrostatic", 0.95),
        DicerPocketResidue(1644, "D", "RNaseIIIb", [42], "electrostatic", 0.9),
        DicerPocketResidue(1648, "E", "RNaseIIIb", [41, 42, 43], "electrostatic", 0.85),
        DicerPocketResidue(1655, "R", "RNaseIIIb", [39, 40], "electrostatic", 0.6),
        DicerPocketResidue(1660, "K", "RNaseIIIb", [43, 44], "electrostatic", 0.5),
        # Helicase domain — contacts dsRNA near 5' end
        DicerPocketResidue(58, "R", "helicase", [0, 1], "electrostatic", 0.4),
        DicerPocketResidue(62, "K", "helicase", [1, 2], "electrostatic", 0.4),
        DicerPocketResidue(95, "R", "helicase", [2, 3], "electrostatic", 0.5),
        DicerPocketResidue(110, "Q", "helicase", [0, 1], "hbond", 0.3),
        DicerPocketResidue(125, "N", "helicase", [3], "hbond", 0.3),
    ]

    if enzyme == "dcr":
        # Fly DCR-1 has slightly different contact residues
        for res in residues:
            res.interaction_strength *= 0.85

    return DicerPocketModel(residues=residues, enzyme=enzyme)
