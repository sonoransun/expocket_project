"""Tests for viewer.genome_editing.off_target."""

from __future__ import annotations

from viewer.genome_editing.off_target import (
    OffTargetScorer,
    compute_specificity,
    find_off_targets_in_dataset,
)
from viewer.genome_editing.sequence_utils import find_pam_sites, get_dual_strand


def test_find_off_targets_returns_list(enriched_dataset):
    # Get a real spacer from the first variant so we have a plausible guide
    variant = enriched_dataset.variants[0]
    sense_dna, _ = get_dual_strand(variant.pre_mirna_sequence)
    pam_sites = find_pam_sites(sense_dna, "NGG", "3prime", 20)
    assert len(pam_sites) > 0, "Need at least one PAM site for this test"
    spacer = pam_sites[0].spacer_seq

    results = find_off_targets_in_dataset(
        guide_spacer=spacer,
        pam_pattern="NGG",
        pam_position="3prime",
        spacer_len=20,
        dataset=enriched_dataset,
        exclude_variant_id=variant.variant,
        threshold=12.0,  # generous threshold to get hits
    )
    assert isinstance(results, list)
    # With 4 variants and a generous threshold, we may find off-targets
    # (not guaranteed, but the function should return a well-formed list)
    for r in results:
        assert "variant_id" in r
        assert "weighted_score" in r
        assert "mismatches" in r


def test_specificity_no_offtargets():
    score = compute_specificity([])
    assert score == 1.0


def test_scorer_facade(enriched_dataset):
    # Get a real spacer from dataset
    variant = enriched_dataset.variants[0]
    sense_dna, _ = get_dual_strand(variant.pre_mirna_sequence)
    pam_sites = find_pam_sites(sense_dna, "NGG", "3prime", 20)
    assert len(pam_sites) > 0
    spacer = pam_sites[0].spacer_seq

    scorer = OffTargetScorer()
    specificity, off_targets = scorer.score(
        guide_spacer=spacer,
        pam_pattern="NGG",
        pam_position="3prime",
        spacer_len=20,
        dataset=enriched_dataset,
        exclude_variant_id=variant.variant,
    )
    assert isinstance(specificity, float)
    assert 0.0 <= specificity <= 1.0
    assert isinstance(off_targets, list)
