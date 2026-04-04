"""Tests for viewer/encoding/nucleotide_properties.py — property vectors and encoding."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.encoding.nucleotide_properties import (
    CANONICAL_PROPERTIES,
    NucleotidePropertyVector,
    encode_sequence,
    get_property,
)


# ---------------------------------------------------------------------------
# CANONICAL_PROPERTIES
# ---------------------------------------------------------------------------

def test_canonical_properties_all_five() -> None:
    """All five canonical nucleotides (A, U, G, C, T) are present."""
    for nt in ("A", "U", "G", "C", "T"):
        assert nt in CANONICAL_PROPERTIES, f"Missing {nt}"


def test_vector_size_is_12() -> None:
    """Static method reports 12 properties."""
    assert NucleotidePropertyVector.vector_size() == 12


def test_as_vector_shape() -> None:
    """as_vector() returns a (12,) float64 array."""
    vec = CANONICAL_PROPERTIES["A"].as_vector()
    assert vec.shape == (12,)
    assert vec.dtype == np.float64


def test_as_vector_values_adenosine() -> None:
    """Spot-check Adenosine MW and vdW volume."""
    a = CANONICAL_PROPERTIES["A"]
    assert a.molecular_weight == pytest.approx(267.24)
    assert a.vdw_volume == pytest.approx(135.0)

    vec = a.as_vector()
    assert vec[0] == pytest.approx(267.24)  # molecular_weight
    assert vec[1] == pytest.approx(135.0)   # vdw_volume


def test_purine_flag() -> None:
    """A and G are purines; U, C, T are not."""
    assert CANONICAL_PROPERTIES["A"].is_purine is True
    assert CANONICAL_PROPERTIES["G"].is_purine is True
    assert CANONICAL_PROPERTIES["U"].is_purine is False
    assert CANONICAL_PROPERTIES["C"].is_purine is False
    assert CANONICAL_PROPERTIES["T"].is_purine is False


# ---------------------------------------------------------------------------
# get_property()
# ---------------------------------------------------------------------------

def test_get_property_canonical() -> None:
    """get_property('A') returns the correct vector object."""
    prop = get_property("A")
    assert prop.symbol == "A"
    assert prop is CANONICAL_PROPERTIES["A"]


def test_get_property_case_insensitive() -> None:
    """get_property('a') works (case-insensitive lookup)."""
    prop = get_property("a")
    assert prop.symbol == "A"
    assert prop.molecular_weight == pytest.approx(267.24)


def test_get_property_unknown_fallback() -> None:
    """get_property('X') falls back to U properties."""
    prop = get_property("X")
    u = CANONICAL_PROPERTIES["U"]
    assert prop.symbol == u.symbol
    assert prop.molecular_weight == pytest.approx(u.molecular_weight)


# ---------------------------------------------------------------------------
# encode_sequence()
# ---------------------------------------------------------------------------

def test_encode_sequence_shape() -> None:
    """encode_sequence('AUGC') returns (4, 12)."""
    mat = encode_sequence("AUGC")
    assert mat.shape == (4, 12)
    assert mat.dtype == np.float64


def test_encode_sequence_empty() -> None:
    """encode_sequence('') returns a zero-length array."""
    mat = encode_sequence("")
    assert mat.shape[0] == 0


def test_encode_sequence_single() -> None:
    """encode_sequence('A') returns (1, 12) matching A vector."""
    mat = encode_sequence("A")
    assert mat.shape == (1, 12)
    expected = CANONICAL_PROPERTIES["A"].as_vector()
    np.testing.assert_array_equal(mat[0], expected)


def test_T_differs_from_U() -> None:
    """T and U have different molecular_weight values."""
    t_mw = CANONICAL_PROPERTIES["T"].molecular_weight
    u_mw = CANONICAL_PROPERTIES["U"].molecular_weight
    assert t_mw != pytest.approx(u_mw), "T and U should have different MW"
