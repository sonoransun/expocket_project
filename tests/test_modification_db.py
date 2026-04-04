"""Tests for viewer/encoding/modification_db.py — chemical modification database."""

from __future__ import annotations

import pytest

from viewer.encoding.modification_db import (
    MODIFICATIONS_DB,
    ChemicalModification,
    applicable_modifications,
    get_modification,
)
from viewer.encoding.nucleotide_properties import PROPERTY_NAMES


# ---------------------------------------------------------------------------
# MODIFICATIONS_DB completeness
# ---------------------------------------------------------------------------

ALL_NINE = ("2OMe", "LNA", "PSI", "m6A", "2F", "PS", "s4U", "m5C", "INO")


def test_all_nine_modifications_present() -> None:
    """All 9 expected modifications are in MODIFICATIONS_DB."""
    for code in ALL_NINE:
        assert code in MODIFICATIONS_DB, f"Missing modification: {code}"
    assert len(MODIFICATIONS_DB) == 9


# ---------------------------------------------------------------------------
# get_modification()
# ---------------------------------------------------------------------------

def test_get_modification_existing() -> None:
    """get_modification returns the correct object for '2OMe'."""
    mod = get_modification("2OMe")
    assert mod is not None
    assert isinstance(mod, ChemicalModification)
    assert mod.code == "2OMe"
    assert mod.full_name == "2'-O-Methylation"


def test_get_modification_nonexistent() -> None:
    """get_modification returns None for an unknown code."""
    assert get_modification("FAKE") is None
    assert get_modification("") is None


# ---------------------------------------------------------------------------
# applicable_modifications()
# ---------------------------------------------------------------------------

def test_applicable_modifications_A() -> None:
    """Check which modifications apply to adenosine."""
    mods = applicable_modifications("A")
    codes = {m.code for m in mods}
    # Universal mods that apply to all 4 nts
    for universal in ("2OMe", "LNA", "2F", "PS"):
        assert universal in codes, f"{universal} should apply to A"
    # A-specific
    assert "m6A" in codes
    assert "INO" in codes
    # U-only / C-only should not apply
    assert "PSI" not in codes
    assert "s4U" not in codes
    assert "m5C" not in codes


def test_applicable_modifications_U() -> None:
    """Check which modifications apply to uridine."""
    mods = applicable_modifications("U")
    codes = {m.code for m in mods}
    for universal in ("2OMe", "LNA", "2F", "PS"):
        assert universal in codes
    assert "PSI" in codes
    assert "s4U" in codes
    # A-only
    assert "m6A" not in codes
    assert "INO" not in codes
    # C-only
    assert "m5C" not in codes


def test_applicable_modifications_G() -> None:
    """Check which modifications apply to guanosine."""
    mods = applicable_modifications("G")
    codes = {m.code for m in mods}
    for universal in ("2OMe", "LNA", "2F", "PS"):
        assert universal in codes
    # G has no specific mods beyond universal ones
    assert "PSI" not in codes
    assert "m6A" not in codes
    assert "m5C" not in codes
    assert "s4U" not in codes
    assert "INO" not in codes


def test_applicable_modifications_C() -> None:
    """Check which modifications apply to cytidine."""
    mods = applicable_modifications("C")
    codes = {m.code for m in mods}
    for universal in ("2OMe", "LNA", "2F", "PS"):
        assert universal in codes
    assert "m5C" in codes
    # Not A or U specific
    assert "m6A" not in codes
    assert "PSI" not in codes
    assert "s4U" not in codes
    assert "INO" not in codes


# ---------------------------------------------------------------------------
# Property deltas validation
# ---------------------------------------------------------------------------

def test_property_deltas_valid_keys() -> None:
    """All delta keys in every modification are valid PROPERTY_NAMES."""
    for code, mod in MODIFICATIONS_DB.items():
        for key in mod.property_deltas:
            assert key in PROPERTY_NAMES, (
                f"Modification '{code}' has invalid delta key '{key}'"
            )


# ---------------------------------------------------------------------------
# Specific applies_to constraints
# ---------------------------------------------------------------------------

def test_PSI_only_U() -> None:
    """PSI (Pseudouridine) applies only to U."""
    psi = MODIFICATIONS_DB["PSI"]
    assert psi.applies_to == ("U",)


def test_m6A_only_A() -> None:
    """m6A (N6-Methyladenosine) applies only to A."""
    m6a = MODIFICATIONS_DB["m6A"]
    assert m6a.applies_to == ("A",)
