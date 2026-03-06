"""Phosphoramidite monomer database for oligonucleotide synthesis planning.

Coupling efficiencies and cost factors based on typical solid-phase
RNA synthesis using 2'-O-TBDMS or 2'-O-TOM protection chemistry.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PhosphoramiditeMonomer:
    """Synthesis monomer for a nucleotide+modification combination."""

    code: str
    full_name: str
    coupling_efficiency: float  # 0.97-0.995 typical
    requires_extended_coupling: bool
    deprotection_group: str  # "Ac", "Bz", "dmf", "iBu", "standard"
    cost_factor: float  # relative to standard RNA monomer (1.0)
    incompatible_adjacent: tuple[str, ...] = ()  # mod codes that conflict nearby


def _build_monomer_db() -> dict[tuple[str, str | None], PhosphoramiditeMonomer]:
    """Build the full monomer database."""
    db: dict[tuple[str, str | None], PhosphoramiditeMonomer] = {}

    # Standard RNA phosphoramidites
    _standard = {
        "A": ("rA-CE", "2'-O-TBDMS Adenosine CE-Phosphoramidite", "Bz"),
        "U": ("rU-CE", "2'-O-TBDMS Uridine CE-Phosphoramidite", "standard"),
        "G": ("rG-CE", "2'-O-TBDMS Guanosine CE-Phosphoramidite", "dmf"),
        "C": ("rC-CE", "2'-O-TBDMS Cytidine CE-Phosphoramidite", "Ac"),
    }
    for nt, (code, name, deprot) in _standard.items():
        db[(nt, None)] = PhosphoramiditeMonomer(
            code=code, full_name=name,
            coupling_efficiency=0.995, requires_extended_coupling=False,
            deprotection_group=deprot, cost_factor=1.0,
        )

    # 2'-O-Methyl (2OMe)
    for nt in ("A", "U", "G", "C"):
        db[(nt, "2OMe")] = PhosphoramiditeMonomer(
            code=f"2OMe-r{nt}-CE", full_name=f"2'-O-Methyl {nt} CE-Phosphoramidite",
            coupling_efficiency=0.990, requires_extended_coupling=False,
            deprotection_group="standard", cost_factor=1.8,
        )

    # Locked Nucleic Acid (LNA)
    for nt in ("A", "U", "G", "C"):
        db[(nt, "LNA")] = PhosphoramiditeMonomer(
            code=f"LNA-{nt}-CE", full_name=f"LNA {nt} CE-Phosphoramidite",
            coupling_efficiency=0.970, requires_extended_coupling=True,
            deprotection_group="standard", cost_factor=3.5,
            incompatible_adjacent=("PS",),
        )

    # 2'-Fluoro (2F)
    for nt in ("A", "U", "G", "C"):
        db[(nt, "2F")] = PhosphoramiditeMonomer(
            code=f"2F-r{nt}-CE", full_name=f"2'-Fluoro {nt} CE-Phosphoramidite",
            coupling_efficiency=0.992, requires_extended_coupling=False,
            deprotection_group="standard", cost_factor=2.0,
        )

    # Phosphorothioate (PS) — backbone modification, applied during oxidation
    for nt in ("A", "U", "G", "C"):
        db[(nt, "PS")] = PhosphoramiditeMonomer(
            code=f"r{nt}-CE+S", full_name=f"{nt} CE-Phosphoramidite (sulfurization)",
            coupling_efficiency=0.980, requires_extended_coupling=False,
            deprotection_group="standard", cost_factor=2.0,
            incompatible_adjacent=("LNA",),
        )

    # Pseudouridine (PSI) — U only
    db[("U", "PSI")] = PhosphoramiditeMonomer(
        code="PSI-CE", full_name="Pseudouridine CE-Phosphoramidite",
        coupling_efficiency=0.988, requires_extended_coupling=False,
        deprotection_group="standard", cost_factor=2.5,
    )

    # N6-Methyladenosine (m6A) — A only
    db[("A", "m6A")] = PhosphoramiditeMonomer(
        code="m6A-CE", full_name="N6-Methyladenosine CE-Phosphoramidite",
        coupling_efficiency=0.985, requires_extended_coupling=False,
        deprotection_group="Bz", cost_factor=3.0,
    )

    # 4-Thiouridine (s4U) — U only
    db[("U", "s4U")] = PhosphoramiditeMonomer(
        code="s4U-CE", full_name="4-Thiouridine CE-Phosphoramidite",
        coupling_efficiency=0.982, requires_extended_coupling=False,
        deprotection_group="standard", cost_factor=3.0,
    )

    # 5-Methylcytosine (m5C) — C only
    db[("C", "m5C")] = PhosphoramiditeMonomer(
        code="m5C-CE", full_name="5-Methylcytosine CE-Phosphoramidite",
        coupling_efficiency=0.990, requires_extended_coupling=False,
        deprotection_group="Ac", cost_factor=2.0,
    )

    # Inosine (INO) — A only (deaminated adenosine)
    db[("A", "INO")] = PhosphoramiditeMonomer(
        code="rI-CE", full_name="Inosine CE-Phosphoramidite",
        coupling_efficiency=0.988, requires_extended_coupling=False,
        deprotection_group="dmf", cost_factor=2.5,
    )

    return db


MONOMER_DB: dict[tuple[str, str | None], PhosphoramiditeMonomer] = _build_monomer_db()


def get_monomer(nucleotide: str, mod_code: str | None = None) -> PhosphoramiditeMonomer:
    """Look up a monomer by nucleotide and optional modification code."""
    nt = nucleotide.upper()
    monomer = MONOMER_DB.get((nt, mod_code))
    if monomer is None and mod_code is not None:
        # Fall back to standard monomer if mod not applicable to this nt
        monomer = MONOMER_DB.get((nt, None))
    if monomer is None:
        # Ultimate fallback: generic U monomer
        monomer = MONOMER_DB[("U", None)]
    return monomer


def check_adjacent_compatibility(mod1: str | None, mod2: str | None) -> bool:
    """Check if two modifications can be adjacent in synthesis.

    Returns True if compatible, False if there's a known conflict.
    """
    if mod1 is None or mod2 is None:
        return True

    # Check both directions — look up any monomer with these mods
    for nt in ("A", "U", "G", "C"):
        m1 = MONOMER_DB.get((nt, mod1))
        if m1 and mod2 in m1.incompatible_adjacent:
            return False
        m2 = MONOMER_DB.get((nt, mod2))
        if m2 and mod1 in m2.incompatible_adjacent:
            return False

    return True
