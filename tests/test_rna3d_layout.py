"""Tests for viewer/rna3d/layout.py — dot-bracket parsing and 3D layout."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.config import PREMIR324_REFERENCE_SEQ, PREMIR324_REFERENCE_STRUCT
from viewer.rna3d.layout import HairpinLayout, compute_hairpin_3d, parse_dot_bracket


# ---------------------------------------------------------------------------
# parse_dot_bracket
# ---------------------------------------------------------------------------

def test_parse_dot_bracket_pairs() -> None:
    """Reference structure produces correct pair indices."""
    pairs, labels = parse_dot_bracket(PREMIR324_REFERENCE_STRUCT)

    # Every pair (i, j) must satisfy: structure[i] == '(' and structure[j] == ')'
    for i, j in pairs:
        assert PREMIR324_REFERENCE_STRUCT[i] == "(", f"position {i} should be '('"
        assert PREMIR324_REFERENCE_STRUCT[j] == ")", f"position {j} should be ')'"

    # Number of pairs must equal the count of '(' characters
    expected_num_pairs = PREMIR324_REFERENCE_STRUCT.count("(")
    assert len(pairs) == expected_num_pairs

    # Each pair index must appear exactly once across all pairs
    all_indices = [idx for pair in pairs for idx in pair]
    assert len(all_indices) == len(set(all_indices)), "pair indices must be unique"


def test_parse_dot_bracket_all_dots() -> None:
    """All dots produce no pairs."""
    pairs, labels = parse_dot_bracket("........")
    assert pairs == []
    assert len(labels) == 8


def test_parse_dot_bracket_empty() -> None:
    """Empty string produces empty result."""
    pairs, labels = parse_dot_bracket("")
    assert pairs == []
    assert labels == []


# ---------------------------------------------------------------------------
# compute_hairpin_3d — shapes
# ---------------------------------------------------------------------------

def test_compute_hairpin_3d_positions_shape() -> None:
    """base_positions has shape (N, 3) for an N-length sequence."""
    seq = PREMIR324_REFERENCE_SEQ
    struct = PREMIR324_REFERENCE_STRUCT
    layout = compute_hairpin_3d(seq, struct)
    n = len(seq)
    assert layout.base_positions.shape == (n, 3)
    assert layout.sequence_length == n


def test_compute_hairpin_3d_cleavage_positions() -> None:
    """DC20-23 positions exist in the cleavage_site_positions dict."""
    seq = PREMIR324_REFERENCE_SEQ
    struct = PREMIR324_REFERENCE_STRUCT
    layout = compute_hairpin_3d(seq, struct)

    for site in [20, 21, 22, 23]:
        assert site in layout.cleavage_site_positions, f"DC{site} missing"
        pos = layout.cleavage_site_positions[site]
        assert pos.shape == (3,), f"DC{site} position should be a 3-vector"


def test_compute_hairpin_3d_backbone_shape() -> None:
    """Backbone spline has more points than the raw sequence (cubic interpolation)."""
    seq = PREMIR324_REFERENCE_SEQ
    struct = PREMIR324_REFERENCE_STRUCT
    layout = compute_hairpin_3d(seq, struct)
    n = len(seq)

    # backbone_spline should be (M, 3) with M > N (since n >= 4, spline uses n*5 points)
    assert layout.backbone_spline.ndim == 2
    assert layout.backbone_spline.shape[1] == 3
    assert layout.backbone_spline.shape[0] > n


def test_compute_hairpin_3d_short_sequence() -> None:
    """A very short sequence (3-4 nt) does not crash."""
    # 3 nt — below the CubicSpline threshold (requires n >= 4)
    layout3 = compute_hairpin_3d("AUG", "(.)")
    assert layout3.base_positions.shape == (3, 3)
    assert layout3.sequence_length == 3
    # backbone falls back to raw positions for n < 4
    assert layout3.backbone_spline.shape[0] == 3

    # 4 nt — exactly at the spline threshold
    layout4 = compute_hairpin_3d("AUGC", "(..)")
    assert layout4.base_positions.shape == (4, 3)
    assert layout4.sequence_length == 4


# ---------------------------------------------------------------------------
# Region labels
# ---------------------------------------------------------------------------

def test_region_labels_count() -> None:
    """Number of region labels equals the sequence length."""
    seq = PREMIR324_REFERENCE_SEQ
    struct = PREMIR324_REFERENCE_STRUCT
    layout = compute_hairpin_3d(seq, struct)
    assert len(layout.region_labels) == len(seq)

    # Each label should be one of the known categories
    valid_labels = {"stem", "loop", "overhang_5p", "overhang_3p", "bulge"}
    for label in layout.region_labels:
        assert label in valid_labels, f"unexpected region label: {label!r}"
