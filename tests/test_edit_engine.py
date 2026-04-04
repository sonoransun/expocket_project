"""Tests for viewer.genome_editing.edit_engine."""

from __future__ import annotations

from viewer.genome_editing.edit_engine import (
    EditEngine,
    apply_deletion,
    apply_insertion,
    apply_snp,
    generate_hdr_template,
    predict_nhej_outcomes,
)


def test_apply_snp():
    assert apply_snp("ACGT", 1, "T") == "ATGT"


def test_apply_snp_boundary():
    # Position 0
    assert apply_snp("ACGT", 0, "G") == "GCGT"
    # Last position
    assert apply_snp("ACGT", 3, "A") == "ACGA"


def test_apply_insertion():
    assert apply_insertion("ACGT", 2, "TT") == "ACTTGT"


def test_apply_deletion():
    assert apply_deletion("ACGT", 1, 2) == "AT"


def test_apply_deletion_past_end():
    # Deleting more than available should not crash; just delete to end
    result = apply_deletion("ACGT", 2, 100)
    assert result == "AC"


def test_hdr_template_contains_edit():
    sense_dna = "A" * 60 + "ACGT" + "A" * 60
    # Edit at position 62 (the 'G' in ACGT), change to 'T'
    template = generate_hdr_template(sense_dna, cut_pos=62, edit_position=62,
                                     new_base="T", flank=50)
    assert "T" in template
    assert len(template) > 0
    assert len(template) <= 100


def test_predict_nhej_outcomes():
    dna = "ACGTACGTACGTACGT"
    outcomes = predict_nhej_outcomes(dna, cut_pos=8, n=5)
    assert isinstance(outcomes, list)
    assert len(outcomes) == 5
    # First outcome is 1-bp deletion at cut position
    assert len(outcomes[0]) == len(dna) - 1
    # Second outcome is 1-bp insertion (duplication)
    assert len(outcomes[1]) == len(dna) + 1


def test_edit_engine_facade_snp():
    engine = EditEngine()
    original_rna = "AUGCAUGCAUGCAUGCAUGC"
    modified_rna, hdr = engine.apply_edit(
        original_rna,
        edit_type="snp",
        edit_position=4,
        edit_sequence="G",
    )
    assert modified_rna != original_rna
    assert modified_rna[4] != original_rna[4]
    assert len(hdr) > 0


def test_edit_engine_facade_unknown_type():
    engine = EditEngine()
    original_rna = "AUGCAUGCAUGCAUGCAUGC"
    modified_rna, hdr = engine.apply_edit(
        original_rna,
        edit_type="unknown_edit_type",
        edit_position=4,
        edit_sequence="G",
    )
    assert modified_rna == original_rna
    assert hdr == ""
