"""Tests for viewer.chemistry.cleavage_predictor.CleavageSitePredictor."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.chemistry.cleavage_predictor import CleavageSitePredictor


class TestPredictorWeights:
    """Tests for predictor construction and weight properties."""

    def test_predictor_weights_exist(self, mock_dataset_256):
        """After construction, _weights has keys 20, 21, 22, 23."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        assert set(predictor._weights.keys()) == {20, 21, 22, 23}

    def test_predictor_weights_shape(self, mock_dataset_256):
        """Each weight vector has shape (48,) matching the 48 summary features."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        for site in [20, 21, 22, 23]:
            assert predictor._weights[site].shape == (48,)


class TestPredictShift:
    """Tests for predict_shift."""

    def test_predict_shift_no_mod(self, mock_dataset_256):
        """Empty modifications dict returns shifts all approximately zero."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        vid = mock_dataset_256.variants[0].variant
        shifts = predictor.predict_shift(vid, modifications={})
        for site in [20, 21, 22, 23]:
            assert abs(shifts[site]) < 1e-6

    def test_predict_shift_with_mod(self, mock_dataset_256):
        """A real modification produces at least one non-zero shift."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        vid = mock_dataset_256.variants[0].variant
        shifts = predictor.predict_shift(vid, modifications={10: "2OMe"})
        any_nonzero = any(abs(shifts[s]) > 1e-10 for s in [20, 21, 22, 23])
        assert any_nonzero

    def test_predict_shift_unknown_variant(self, mock_dataset_256):
        """Unknown variant returns all-zero shifts."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        shifts = predictor.predict_shift("NONEXISTENT", modifications={5: "2OMe"})
        for site in [20, 21, 22, 23]:
            assert shifts[site] == 0.0


class TestPredictDcRatioShift:
    """Tests for predict_dc_ratio_shift."""

    def test_predict_dc_ratio_shift(self, mock_dataset_256):
        """dc_ratio_shift equals shifts[21] - shifts[22]."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        vid = mock_dataset_256.variants[0].variant
        mods = {10: "2OMe"}
        shifts = predictor.predict_shift(vid, mods)
        ratio_shift = predictor.predict_dc_ratio_shift(vid, mods)
        assert abs(ratio_shift - (shifts[21] - shifts[22])) < 1e-10


class TestPredictBaseChange:
    """Tests for predict_base_change."""

    def test_predict_base_change_valid(self, mock_dataset_256):
        """A nucleotide substitution produces non-zero shifts."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        vid = mock_dataset_256.variants[0].variant
        # Substitute position 5 (a C in the reference) with G
        shifts = predictor.predict_base_change(vid, position=5, new_base="G")
        any_nonzero = any(abs(shifts[s]) > 1e-10 for s in [20, 21, 22, 23])
        assert any_nonzero


class TestPredictShiftBatch:
    """Tests for predict_shift_batch."""

    def test_predict_shift_batch_empty(self, mock_dataset_256):
        """Empty modifications_list returns empty list."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        vid = mock_dataset_256.variants[0].variant
        result = predictor.predict_shift_batch(vid, [])
        assert result == []

    def test_predict_shift_batch_single(self, mock_dataset_256):
        """Single-element batch matches predict_shift result."""
        predictor = CleavageSitePredictor(mock_dataset_256)
        vid = mock_dataset_256.variants[0].variant
        mods = {10: "2OMe"}
        single = predictor.predict_shift(vid, mods)
        batch = predictor.predict_shift_batch(vid, [mods])
        assert len(batch) == 1
        for site in [20, 21, 22, 23]:
            assert abs(batch[0][site] - single[site]) < 1e-10


class TestSummarizeSingle:
    """Tests for _summarize_single static method."""

    def test_summarize_single_shape(self):
        """Returns a 48-element feature vector."""
        props = np.random.rand(40, 12)
        result = CleavageSitePredictor._summarize_single(props)
        assert result.shape == (48,)
