"""Tests for viewer/data/schema.py — dataclasses and dataset operations."""

from __future__ import annotations

import logging

import numpy as np
import pytest

from viewer.data.schema import (
    CleavageRecord,
    EnrichedVariantDataset,
    ModificationState,
    ReplacementResult,
    SynthesisPlan,
    VariantDataset,
    VariantInfo,
)


# ---------------------------------------------------------------------------
# VariantInfo
# ---------------------------------------------------------------------------

def test_variant_info_construction(sample_variant_info: VariantInfo) -> None:
    """Fields are stored correctly from fixture."""
    v = sample_variant_info
    assert v.variant == "A_AAA"
    assert v.group == "A"
    assert v.randomized_nts == "AAA"
    assert len(v.pre_mirna_sequence) > 0
    assert v.flanking_length_5p == 0
    assert v.concrete_struct != ""


def test_variant_info_empty_sequence_allowed(caplog: pytest.LogCaptureFixture) -> None:
    """An empty sequence is accepted but triggers a warning log."""
    with caplog.at_level(logging.WARNING, logger="viewer.data.schema"):
        v = VariantInfo(
            variant="test",
            group="A",
            randomized_nts="AAA",
            pre_mirna_sequence="",
            concrete_struct=".",
            flanking_length_5p=0,
        )
    assert v.pre_mirna_sequence == ""
    assert any("empty sequence" in msg for msg in caplog.messages)


def test_variant_info_negative_flanking_clamped(caplog: pytest.LogCaptureFixture) -> None:
    """Negative flanking_length_5p is clamped to 0 with a warning."""
    with caplog.at_level(logging.WARNING, logger="viewer.data.schema"):
        v = VariantInfo(
            variant="neg",
            group="G",
            randomized_nts="GGG",
            pre_mirna_sequence="AUGC",
            concrete_struct="(..)",
            flanking_length_5p=-5,
        )
    assert v.flanking_length_5p == 0
    assert any("negative flanking_length_5p" in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# CleavageRecord
# ---------------------------------------------------------------------------

def test_cleavage_record_defaults() -> None:
    """Optional float fields default to 0.0."""
    rec = CleavageRecord(variant="v1", cleavage_site=21, mean_accuracy=0.5)
    assert rec.mean_positional_efficiency == 0.0
    assert rec.mean_global_efficiency == 0.0
    assert rec.accuracy_rep1 == 0.0
    assert rec.accuracy_rep2 == 0.0
    assert rec.accuracy_rep3 == 0.0


def test_cleavage_record_unusual_site_allowed(caplog: pytest.LogCaptureFixture) -> None:
    """A non-standard site (999) is accepted but logged as unusual."""
    with caplog.at_level(logging.WARNING, logger="viewer.data.schema"):
        rec = CleavageRecord(variant="v1", cleavage_site=999, mean_accuracy=0.3)
    assert rec.cleavage_site == 999
    assert any("unusual cleavage_site=999" in msg for msg in caplog.messages)


def test_cleavage_record_accuracy_clamped() -> None:
    """mean_accuracy > 1.0 is clamped to 1.0 and < 0.0 to 0.0."""
    high = CleavageRecord(variant="v1", cleavage_site=21, mean_accuracy=1.5)
    assert high.mean_accuracy == 1.0

    low = CleavageRecord(variant="v1", cleavage_site=21, mean_accuracy=-0.3)
    assert low.mean_accuracy == 0.0


# ---------------------------------------------------------------------------
# VariantDataset
# ---------------------------------------------------------------------------

def test_variant_dataset_variant_ids(small_dataset: VariantDataset) -> None:
    """variant_ids returns ordered list of variant names."""
    ids = small_dataset.variant_ids
    assert ids == ["A_AAA", "T_TTT", "G_GGG", "C_CCC"]


def test_variant_dataset_get_variant_found(small_dataset: VariantDataset) -> None:
    """get_variant returns VariantInfo for a known variant."""
    v = small_dataset.get_variant("G_GGG")
    assert v is not None
    assert v.variant == "G_GGG"
    assert v.group == "G"


def test_variant_dataset_get_variant_missing(small_dataset: VariantDataset) -> None:
    """get_variant returns None for an unknown variant."""
    assert small_dataset.get_variant("NONEXISTENT") is None


def test_variant_dataset_get_cleavage_found(small_dataset: VariantDataset) -> None:
    """get_cleavage returns CleavageRecord for a known variant and site."""
    rec = small_dataset.get_cleavage("A_AAA", 21)
    assert rec is not None
    assert rec.cleavage_site == 21
    assert rec.mean_accuracy == pytest.approx(0.55)


def test_variant_dataset_get_cleavage_missing(small_dataset: VariantDataset) -> None:
    """get_cleavage returns None for unknown variant or site."""
    assert small_dataset.get_cleavage("MISSING", 21) is None
    assert small_dataset.get_cleavage("A_AAA", 99) is None


def test_variant_dataset_empty(empty_dataset: VariantDataset) -> None:
    """Empty dataset has empty variant_ids."""
    assert empty_dataset.variant_ids == []
    assert empty_dataset.get_variant("x") is None
    assert empty_dataset.get_cleavage("x", 21) is None


# ---------------------------------------------------------------------------
# ModificationState
# ---------------------------------------------------------------------------

def test_modification_state_apply_remove_clear() -> None:
    """Full CRUD cycle on ModificationState."""
    ms = ModificationState(variant="v1")
    assert ms.modifications == {}

    # Apply
    ms.apply(5, "2OMe")
    ms.apply(10, "LNA")
    assert ms.modifications == {5: "2OMe", 10: "LNA"}

    # Update (overwrite)
    ms.apply(5, "PSI")
    assert ms.modifications[5] == "PSI"

    # Remove
    ms.remove(10)
    assert 10 not in ms.modifications

    # Remove non-existent is safe
    ms.remove(99)

    # Clear
    ms.clear()
    assert ms.modifications == {}


# ---------------------------------------------------------------------------
# EnrichedVariantDataset
# ---------------------------------------------------------------------------

def test_enriched_dataset_get_modification_state(
    enriched_dataset: EnrichedVariantDataset,
) -> None:
    """get_modification_state auto-creates on first access."""
    ms = enriched_dataset.get_modification_state("A_AAA")
    assert isinstance(ms, ModificationState)
    assert ms.variant == "A_AAA"
    assert ms.modifications == {}

    # Same object returned on second call
    ms2 = enriched_dataset.get_modification_state("A_AAA")
    assert ms is ms2

    # Auto-creates for novel variant id
    ms3 = enriched_dataset.get_modification_state("NOVEL_ID")
    assert ms3.variant == "NOVEL_ID"


# ---------------------------------------------------------------------------
# SynthesisPlan and ReplacementResult defaults
# ---------------------------------------------------------------------------

def test_synthesis_plan_defaults() -> None:
    """SynthesisPlan fields have expected defaults."""
    sp = SynthesisPlan(variant_id="v1", sequence="AUGC")
    assert sp.modifications == {}
    assert sp.steps == []
    assert sp.total_yield == 0.0
    assert sp.total_cost_factor == 0.0
    assert sp.incompatibilities == []
    assert sp.scale_nmol == 200


def test_replacement_result_defaults() -> None:
    """ReplacementResult fields have expected defaults."""
    rr = ReplacementResult()
    assert rr.positions == []
    assert rr.original_bases == []
    assert rr.replacements == []
    assert rr.shifts == {}
    assert rr.dc_ratio_shift == 0.0
    assert rr.synthesis_feasible is True
    assert rr.estimated_yield == 0.0
