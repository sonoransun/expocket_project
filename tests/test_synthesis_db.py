"""Tests for viewer.encoding.synthesis_db."""

from __future__ import annotations

import pytest

from viewer.encoding.synthesis_db import (
    MONOMER_DB,
    check_adjacent_compatibility,
    get_monomer,
)


# ---------------------------------------------------------------------------
# Standard monomers
# ---------------------------------------------------------------------------

def test_standard_monomers_all_four():
    for nt in ("A", "U", "G", "C"):
        assert (nt, None) in MONOMER_DB


def test_standard_coupling_efficiency():
    for nt in ("A", "U", "G", "C"):
        monomer = MONOMER_DB[(nt, None)]
        assert monomer.coupling_efficiency == 0.995


# ---------------------------------------------------------------------------
# get_monomer
# ---------------------------------------------------------------------------

def test_get_monomer_standard():
    monomer = get_monomer("A")
    assert monomer.code == "rA-CE"
    assert monomer.coupling_efficiency == 0.995


def test_get_monomer_with_mod():
    monomer = get_monomer("A", "2OMe")
    assert "2OMe" in monomer.code
    assert monomer.coupling_efficiency < 0.995


def test_get_monomer_inapplicable_fallback():
    # PSI only applies to U, so get_monomer("G", "PSI") should fall back to standard G
    monomer = get_monomer("G", "PSI")
    assert monomer.code == "rG-CE"
    assert monomer.coupling_efficiency == 0.995


def test_get_monomer_unknown_nt():
    # "X" is not a known nucleotide; should fall back to standard U
    monomer = get_monomer("X")
    assert monomer.code == "rU-CE"


# ---------------------------------------------------------------------------
# check_adjacent_compatibility
# ---------------------------------------------------------------------------

def test_check_compatibility_LNA_PS():
    # LNA and PS are declared incompatible in both directions
    assert check_adjacent_compatibility("LNA", "PS") is False
    assert check_adjacent_compatibility("PS", "LNA") is False


def test_check_compatibility_both_none():
    assert check_adjacent_compatibility(None, None) is True


# ---------------------------------------------------------------------------
# Modified monomer properties
# ---------------------------------------------------------------------------

def test_modified_lower_efficiency():
    for key, monomer in MONOMER_DB.items():
        nt, mod = key
        if mod is not None:
            assert monomer.coupling_efficiency < 0.995, (
                f"Modified monomer {key} should have coupling_efficiency < 0.995"
            )


def test_LNA_extended_coupling():
    for nt in ("A", "U", "G", "C"):
        monomer = MONOMER_DB[(nt, "LNA")]
        assert monomer.requires_extended_coupling is True
