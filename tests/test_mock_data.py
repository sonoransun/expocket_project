"""Tests for viewer/data/mock_data.py — mock dataset generation and enrichment."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.config import GROUPS
from viewer.data.mock_data import enrich_dataset, generate_mock_dataset
from viewer.data.schema import EnrichedVariantDataset, VariantDataset


# ---------------------------------------------------------------------------
# generate_mock_dataset — count and grouping
# ---------------------------------------------------------------------------

def test_generate_mock_dataset_256() -> None:
    """Default mock dataset returns exactly 256 variants."""
    ds = generate_mock_dataset(seed=42)
    assert isinstance(ds, VariantDataset)
    assert len(ds.variants) == 256


def test_generate_mock_dataset_groups() -> None:
    """Each of the 4 groups (A, T, G, C) contains exactly 64 variants."""
    ds = generate_mock_dataset(seed=42)
    for group in GROUPS:
        count = sum(1 for v in ds.variants if v.group == group)
        assert count == 64, f"group {group} has {count} variants, expected 64"


# ---------------------------------------------------------------------------
# Determinism and seed control
# ---------------------------------------------------------------------------

def test_generate_mock_dataset_deterministic() -> None:
    """Same seed produces identical datasets."""
    ds1 = generate_mock_dataset(seed=99)
    ds2 = generate_mock_dataset(seed=99)

    assert len(ds1.variants) == len(ds2.variants)
    for v1, v2 in zip(ds1.variants, ds2.variants):
        assert v1.variant == v2.variant
        assert v1.pre_mirna_sequence == v2.pre_mirna_sequence

    # Cleavage accuracies must match exactly
    for vid in ds1.variant_ids:
        for site in [20, 21, 22, 23]:
            r1 = ds1.get_cleavage(vid, site)
            r2 = ds2.get_cleavage(vid, site)
            assert r1 is not None and r2 is not None
            assert r1.mean_accuracy == pytest.approx(r2.mean_accuracy)


def test_generate_mock_dataset_different_seeds() -> None:
    """Different seeds produce different cleavage data."""
    ds1 = generate_mock_dataset(seed=1)
    ds2 = generate_mock_dataset(seed=2)

    # At least some accuracies should differ
    diffs = 0
    for vid in ds1.variant_ids:
        r1 = ds1.get_cleavage(vid, 21)
        r2 = ds2.get_cleavage(vid, 21)
        if r1 is not None and r2 is not None:
            if abs(r1.mean_accuracy - r2.mean_accuracy) > 1e-6:
                diffs += 1
    assert diffs > 0, "different seeds should produce different data"


# ---------------------------------------------------------------------------
# Cleavage accuracy properties
# ---------------------------------------------------------------------------

def test_cleavage_accuracies_sum_to_one() -> None:
    """DC20+DC21+DC22+DC23 accuracies sum to approximately 1.0 for each variant."""
    ds = generate_mock_dataset(seed=42)
    for vid in ds.variant_ids:
        total = 0.0
        for site in [20, 21, 22, 23]:
            rec = ds.get_cleavage(vid, site)
            assert rec is not None, f"missing cleavage for {vid} site {site}"
            total += rec.mean_accuracy
        assert total == pytest.approx(1.0, abs=0.01), (
            f"variant {vid}: accuracy sum = {total:.4f}, expected ~1.0"
        )


def test_cleavage_accuracies_in_range() -> None:
    """All cleavage accuracies are in [0, 1]."""
    ds = generate_mock_dataset(seed=42)
    for vid in ds.variant_ids:
        for site in [20, 21, 22, 23]:
            rec = ds.get_cleavage(vid, site)
            assert rec is not None
            assert 0.0 <= rec.mean_accuracy <= 1.0, (
                f"variant {vid} site {site}: accuracy {rec.mean_accuracy} out of range"
            )


# ---------------------------------------------------------------------------
# enrich_dataset
# ---------------------------------------------------------------------------

def test_enrich_dataset_adds_property_matrix() -> None:
    """Enriched dataset has a property_matrix with shape (256, L, 12)."""
    ds = generate_mock_dataset(seed=42)
    enriched = enrich_dataset(ds)
    assert isinstance(enriched, EnrichedVariantDataset)
    assert enriched.property_matrix is not None
    pm = enriched.property_matrix
    assert pm.ndim == 3
    assert pm.shape[0] == 256
    assert pm.shape[2] == 12  # 12 physicochemical properties


def test_enrich_dataset_adds_summary_features() -> None:
    """Enriched dataset has summary_features with shape (256, 48)."""
    ds = generate_mock_dataset(seed=42)
    enriched = enrich_dataset(ds)
    assert enriched.summary_features is not None
    assert enriched.summary_features.shape == (256, 48)


def test_enrich_dataset_adds_dicer_pocket() -> None:
    """Enriched dataset has a non-None dicer_pocket."""
    ds = generate_mock_dataset(seed=42)
    enriched = enrich_dataset(ds)
    assert enriched.dicer_pocket is not None


# ---------------------------------------------------------------------------
# Variant sequence diversity
# ---------------------------------------------------------------------------

def test_mock_variant_sequences_differ() -> None:
    """Not all 256 variant sequences are identical."""
    ds = generate_mock_dataset(seed=42)
    sequences = {v.pre_mirna_sequence for v in ds.variants}
    # Each variant has unique 5' nt + 3' 3-nt combo, so sequences should vary
    assert len(sequences) > 1, "all 256 variants have the same sequence"
