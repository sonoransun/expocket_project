"""Tests for viewer.chemistry.modification_engine.ModificationEngine."""

from __future__ import annotations

import numpy as np
import pytest

from viewer.chemistry.modification_engine import ModificationEngine


class TestApplyModification:
    """Tests for ModificationEngine.apply_modification."""

    def test_apply_modification_success(self, enriched_dataset):
        """Applying 2OMe at a valid position succeeds and records the mod."""
        engine = ModificationEngine(enriched_dataset)
        state = engine.apply_modification("A_AAA", 0, "2OMe")
        assert state.modifications[0] == "2OMe"

    def test_apply_modification_unknown_variant(self, enriched_dataset):
        """Applying a modification to an unknown variant raises ValueError."""
        engine = ModificationEngine(enriched_dataset)
        with pytest.raises(ValueError, match="Unknown variant"):
            engine.apply_modification("NONEXISTENT", 0, "2OMe")

    def test_apply_modification_position_out_of_range(self, enriched_dataset):
        """A position beyond the sequence length raises ValueError."""
        engine = ModificationEngine(enriched_dataset)
        with pytest.raises(ValueError, match="out of range"):
            engine.apply_modification("A_AAA", 999, "2OMe")

    def test_apply_modification_unknown_mod(self, enriched_dataset):
        """An unknown modification code raises ValueError."""
        engine = ModificationEngine(enriched_dataset)
        with pytest.raises(ValueError, match="Unknown modification"):
            engine.apply_modification("A_AAA", 0, "FAKE_MOD")

    def test_apply_modification_incompatible_nt(self, enriched_dataset):
        """PSI applied to a G position raises ValueError (PSI only applies to U)."""
        engine = ModificationEngine(enriched_dataset)
        # Position 1 in A_AAA is G (from reference: xGCAUCCCCU...)
        with pytest.raises(ValueError, match="cannot be applied"):
            engine.apply_modification("A_AAA", 1, "PSI")


class TestRemoveAndClear:
    """Tests for remove_modification and clear_modifications."""

    def test_remove_modification(self, enriched_dataset):
        """Apply then remove leaves the position unmodified."""
        engine = ModificationEngine(enriched_dataset)
        engine.apply_modification("A_AAA", 0, "2OMe")
        state = engine.remove_modification("A_AAA", 0)
        assert 0 not in state.modifications

    def test_clear_modifications(self, enriched_dataset):
        """clear_modifications empties all modifications for the variant."""
        engine = ModificationEngine(enriched_dataset)
        engine.apply_modification("A_AAA", 0, "2OMe")
        engine.apply_modification("A_AAA", 4, "PSI")  # pos 4 is U
        state = engine.clear_modifications("A_AAA")
        assert len(state.modifications) == 0


class TestGetModifiedProperties:
    """Tests for get_modified_properties."""

    def test_get_modified_properties_shape(self, enriched_dataset):
        """Returns (seq_len, 12) property matrix."""
        engine = ModificationEngine(enriched_dataset)
        variant = enriched_dataset.get_variant("A_AAA")
        seq_len = len(variant.pre_mirna_sequence)
        props = engine.get_modified_properties("A_AAA")
        assert props.shape == (seq_len, 12)


class TestApplicableAtPosition:
    """Tests for applicable_at_position."""

    def test_applicable_at_position(self, enriched_dataset):
        """Position 0 of A_AAA is A; should include 2OMe, LNA, m6A, 2F, PS, INO."""
        engine = ModificationEngine(enriched_dataset)
        codes = engine.applicable_at_position("A_AAA", 0)
        # A can take: 2OMe, LNA, m6A, 2F, PS, INO
        assert "2OMe" in codes
        assert "LNA" in codes
        assert "m6A" in codes
        assert "2F" in codes
        assert "PS" in codes
        assert "INO" in codes
        # PSI and s4U apply only to U, m5C only to C
        assert "PSI" not in codes
        assert "s4U" not in codes
        assert "m5C" not in codes
