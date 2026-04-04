"""Tests for viewer.genome_editing.impact_predictor."""

from __future__ import annotations

from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
from viewer.genome_editing.edit_engine import rna_sequence_after_snp
from viewer.genome_editing.impact_predictor import EditImpactPredictor, _structural_note


def test_predict_no_change(enriched_dataset):
    predictor = CleavageSitePredictor(enriched_dataset)
    impact_pred = EditImpactPredictor(predictor, enriched_dataset)

    variant = enriched_dataset.variants[0]
    # Pass the same RNA as "modified" -- deltas should be ~0
    result = impact_pred.predict(variant.variant, variant.pre_mirna_sequence)
    for site in [20, 21, 22, 23]:
        assert abs(result.delta[site]) < 1e-6


def test_predict_snp_in_cleavage_zone(enriched_dataset):
    predictor = CleavageSitePredictor(enriched_dataset)
    impact_pred = EditImpactPredictor(predictor, enriched_dataset)

    variant = enriched_dataset.variants[0]
    # Apply a SNP at position 20 (within cleavage zone)
    modified_rna = rna_sequence_after_snp(variant.pre_mirna_sequence, 20, "G")
    result = impact_pred.predict(variant.variant, modified_rna)

    # At least one cleavage site should have a non-zero delta
    any_nonzero = any(abs(result.delta[s]) > 1e-9 for s in [20, 21, 22, 23])
    assert any_nonzero


def test_predict_unknown_variant(enriched_dataset):
    predictor = CleavageSitePredictor(enriched_dataset)
    impact_pred = EditImpactPredictor(predictor, enriched_dataset)

    result = impact_pred.predict("NONEXISTENT_999", "AUGCAUGC")
    for site in [20, 21, 22, 23]:
        assert result.delta[site] == 0.0


def test_structural_note_no_change():
    note = _structural_note("AUGCAUGCAUGC", "AUGCAUGCAUGC")
    assert "no change" in note.lower()
