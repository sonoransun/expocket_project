"""Tests for viewer.encoding.property_calculator."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.encoding.property_calculator import (
    SUMMARY_FEATURE_NAMES,
    compute_dataset_property_matrix,
    compute_modified_properties,
    compute_summary_features,
    compute_variant_properties,
)


# ---------------------------------------------------------------------------
# compute_variant_properties
# ---------------------------------------------------------------------------

def test_compute_variant_properties_shape():
    seq = "AUGCAUGCAU"
    result = compute_variant_properties(seq)
    assert result.shape == (len(seq), 12)


# ---------------------------------------------------------------------------
# compute_modified_properties
# ---------------------------------------------------------------------------

def test_compute_modified_properties_no_mods():
    seq = "AUGC"
    unmodified = compute_variant_properties(seq)
    modified = compute_modified_properties(seq, {})
    np.testing.assert_array_equal(modified, unmodified)


def test_compute_modified_properties_single_mod():
    # Sequence starts with A; 2OMe adds +14.0 to molecular_weight (index 0).
    seq = "AUGC"
    modified = compute_modified_properties(seq, {0: "2OMe"})
    unmodified = compute_variant_properties(seq)
    # MW at position 0 should be higher after 2OMe
    assert modified[0, 0] > unmodified[0, 0]
    # Positions 1-3 should be unchanged
    np.testing.assert_array_equal(modified[1:], unmodified[1:])


def test_compute_modified_properties_invalid_pos_ignored():
    seq = "AUGC"
    unmodified = compute_variant_properties(seq)
    # Position -1 and position >= len(seq) should be silently ignored
    modified = compute_modified_properties(seq, {-1: "2OMe", len(seq): "2OMe"})
    np.testing.assert_array_equal(modified, unmodified)


def test_compute_modified_properties_unknown_mod_ignored():
    seq = "AUGC"
    unmodified = compute_variant_properties(seq)
    modified = compute_modified_properties(seq, {0: "FAKE"})
    np.testing.assert_array_equal(modified, unmodified)


# ---------------------------------------------------------------------------
# compute_dataset_property_matrix
# ---------------------------------------------------------------------------

def test_compute_dataset_property_matrix_shape(small_dataset):
    result = compute_dataset_property_matrix(small_dataset)
    max_len = max(len(v.pre_mirna_sequence) for v in small_dataset.variants)
    assert result.shape == (4, max_len, 12)


def test_compute_dataset_property_matrix_empty(empty_dataset):
    result = compute_dataset_property_matrix(empty_dataset)
    assert result.shape == (0, 0, 12)


# ---------------------------------------------------------------------------
# compute_summary_features
# ---------------------------------------------------------------------------

def test_compute_summary_features_shape(small_dataset):
    result = compute_summary_features(small_dataset)
    assert result.shape == (4, 48)


def test_compute_summary_features_empty(empty_dataset):
    result = compute_summary_features(empty_dataset)
    assert result.shape == (0, 48)


# ---------------------------------------------------------------------------
# SUMMARY_FEATURE_NAMES
# ---------------------------------------------------------------------------

def test_summary_feature_names_length():
    assert len(SUMMARY_FEATURE_NAMES) == 48
