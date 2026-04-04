"""Tests for viewer.chemistry.double_screen.DoubleReplacementScreener."""

from __future__ import annotations

import pytest

from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
from viewer.chemistry.double_screen import DoubleReplacementScreener
from viewer.chemistry.modification_engine import ModificationEngine


class TestScreenDouble:
    """Tests for DoubleReplacementScreener.screen_double."""

    def test_screen_double_returns_results(self, mock_dataset_256):
        """Screening a known variant returns a non-empty results list."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        engine = ModificationEngine(mock_dataset_256)
        screener = DoubleReplacementScreener(mock_dataset_256, predictor, engine)
        vid = mock_dataset_256.variants[0].variant
        results = screener.screen_double(vid, top_n=10)
        assert len(results) > 0

    def test_screen_double_unknown_empty(self, mock_dataset_256):
        """Screening an unknown variant returns an empty list."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        engine = ModificationEngine(mock_dataset_256)
        screener = DoubleReplacementScreener(mock_dataset_256, predictor, engine)
        results = screener.screen_double("NONEXISTENT")
        assert results == []

    def test_screen_double_top_n(self, mock_dataset_256):
        """screen_double respects the top_n limit."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        engine = ModificationEngine(mock_dataset_256)
        screener = DoubleReplacementScreener(mock_dataset_256, predictor, engine)
        vid = mock_dataset_256.variants[0].variant
        results = screener.screen_double(vid, top_n=3)
        assert len(results) <= 3


class TestSynergyCalculation:
    """Tests for synergy via compare_single_vs_double."""

    def test_synergy_calculation(self, mock_dataset_256):
        """synergy == delta_ratio - (single1_ratio + single2_ratio)."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        engine = ModificationEngine(mock_dataset_256)
        screener = DoubleReplacementScreener(mock_dataset_256, predictor, engine)
        vid = mock_dataset_256.variants[0].variant

        # Use two distinct positions that both accept 2OMe
        # Position 10 and 20 should both be valid positions in a 63-nt sequence
        comparison = screener.compare_single_vs_double(
            vid, pos1=10, mod1="2OMe", pos2=20, mod2="2OMe"
        )
        single1_ratio = comparison["single_1"][21] - comparison["single_1"][22]
        single2_ratio = comparison["single_2"][21] - comparison["single_2"][22]
        double_ratio = comparison["double"][21] - comparison["double"][22]
        expected_synergy = double_ratio - (single1_ratio + single2_ratio)
        assert abs(comparison["synergy"] - expected_synergy) < 1e-10


class TestEstimateYield:
    """Tests for _estimate_yield static method."""

    def test_estimate_yield(self, mock_dataset_256):
        """Estimated yield is positive and at most 1.0."""
        variant = mock_dataset_256.variants[0]
        seq = variant.pre_mirna_sequence
        yield_val = DoubleReplacementScreener._estimate_yield(seq, {10: "2OMe"})
        assert yield_val > 0.0
        assert yield_val <= 1.0
