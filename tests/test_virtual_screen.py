"""Tests for viewer.chemistry.virtual_screen.VirtualScreener."""

from __future__ import annotations

import pytest

from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
from viewer.chemistry.modification_engine import ModificationEngine
from viewer.chemistry.virtual_screen import VirtualScreener


class TestScreenVariant:
    """Tests for VirtualScreener.screen_variant."""

    def test_screen_returns_results(self, mock_dataset_256):
        """Screening a known variant returns a non-empty results list."""
        engine = ModificationEngine(mock_dataset_256)
        predictor = CleavageSitePredictor(mock_dataset_256)
        screener = VirtualScreener(mock_dataset_256, engine, predictor)
        vid = mock_dataset_256.variants[0].variant
        results = screener.screen_variant(vid)
        assert len(results) > 0

    def test_screen_unknown_empty(self, mock_dataset_256):
        """Screening an unknown variant returns an empty list."""
        engine = ModificationEngine(mock_dataset_256)
        predictor = CleavageSitePredictor(mock_dataset_256)
        screener = VirtualScreener(mock_dataset_256, engine, predictor)
        results = screener.screen_variant("NONEXISTENT")
        assert results == []

    def test_screen_result_fields(self, mock_dataset_256):
        """Each result has delta_ratio == delta_dc21 - delta_dc22."""
        engine = ModificationEngine(mock_dataset_256)
        predictor = CleavageSitePredictor(mock_dataset_256)
        screener = VirtualScreener(mock_dataset_256, engine, predictor)
        vid = mock_dataset_256.variants[0].variant
        results = screener.screen_variant(vid)
        for r in results[:10]:  # check first 10 to keep test fast
            assert abs(r.delta_ratio - (r.delta_dc21 - r.delta_dc22)) < 1e-12


class TestRankTopN:
    """Tests for rank_by_dc_ratio_shift."""

    def test_rank_top_n(self, mock_dataset_256):
        """rank_by_dc_ratio_shift returns at most top_n results."""
        engine = ModificationEngine(mock_dataset_256)
        predictor = CleavageSitePredictor(mock_dataset_256)
        screener = VirtualScreener(mock_dataset_256, engine, predictor)
        vid = mock_dataset_256.variants[0].variant
        results = screener.rank_by_dc_ratio_shift(vid, top_n=5)
        assert len(results) <= 5
