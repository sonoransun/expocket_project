"""Tests for viewer.encoding.protein_descriptors."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.encoding.protein_descriptors import (
    AA_PROPERTIES,
    AminoAcidProperties,
    DicerPocketModel,
    build_mock_dicer_pocket,
)


# ---------------------------------------------------------------------------
# AA_PROPERTIES
# ---------------------------------------------------------------------------

STANDARD_AAS = set("ARNDCEQGHILKMFPSTWYV")


def test_aa_properties_all_20():
    assert set(AA_PROPERTIES.keys()) == STANDARD_AAS


def test_aa_as_vector_shape():
    for aa, props in AA_PROPERTIES.items():
        vec = props.as_vector()
        assert vec.shape == (5,), f"as_vector() shape for {aa} should be (5,)"
        assert vec.dtype == np.float64


# ---------------------------------------------------------------------------
# build_mock_dicer_pocket — human DICER
# ---------------------------------------------------------------------------

def test_dicer_pocket_model_hdicer():
    pocket = build_mock_dicer_pocket("hdicer")
    assert pocket.enzyme == "hdicer"
    assert len(pocket.residues) == 30


def test_dicer_pocket_model_hdicer_domains():
    pocket = build_mock_dicer_pocket("hdicer")
    assert len(pocket.domains) == 5


# ---------------------------------------------------------------------------
# build_mock_dicer_pocket — fly DCR-1
# ---------------------------------------------------------------------------

def test_dicer_pocket_model_dcr():
    pocket_h = build_mock_dicer_pocket("hdicer")
    pocket_d = build_mock_dicer_pocket("dcr")
    assert pocket_d.enzyme == "dcr"
    assert len(pocket_d.residues) == 30
    # Fly DCR-1 has scaled-down interaction strengths (0.85x)
    for rh, rd in zip(pocket_h.residues, pocket_d.residues):
        assert rd.interaction_strength == pytest.approx(
            rh.interaction_strength * 0.85, abs=1e-9
        )


# ---------------------------------------------------------------------------
# contact_matrix
# ---------------------------------------------------------------------------

def test_contact_matrix_shape():
    pocket = build_mock_dicer_pocket("hdicer")
    seq_len = 63
    mat = pocket.contact_matrix(seq_len)
    assert mat.shape == (30, seq_len)


# ---------------------------------------------------------------------------
# domain_residues
# ---------------------------------------------------------------------------

def test_domain_residues():
    pocket = build_mock_dicer_pocket("hdicer")
    expected_counts = {
        "PAZ": 8,
        "platform": 4,
        "RNaseIIIa": 7,
        "RNaseIIIb": 6,
        "helicase": 5,
    }
    for domain, count in expected_counts.items():
        assert len(pocket.domain_residues(domain)) == count, (
            f"Domain {domain} should have {count} residues"
        )


# ---------------------------------------------------------------------------
# domains property
# ---------------------------------------------------------------------------

def test_domains_property():
    pocket = build_mock_dicer_pocket("hdicer")
    domains = pocket.domains
    assert len(domains) == 5
    # Sorted alphabetically
    assert domains == sorted(domains)
