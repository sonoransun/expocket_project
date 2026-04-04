"""Tests for viewer.genome_editing.sequence_utils."""

from __future__ import annotations

from viewer.genome_editing.sequence_utils import (
    count_mismatches_weighted,
    dna_to_rna,
    find_pam_sites,
    gc_content,
    get_dual_strand,
    has_homopolymer,
    reverse_complement_dna,
    reverse_complement_rna,
    rna_to_dna,
)


def test_rna_to_dna():
    assert rna_to_dna("AUGC") == "ATGC"


def test_dna_to_rna():
    assert dna_to_rna("ATGC") == "AUGC"


def test_reverse_complement_dna():
    assert reverse_complement_dna("ATGC") == "GCAT"


def test_reverse_complement_rna():
    assert reverse_complement_rna("AUGC") == "GCAU"


def test_get_dual_strand():
    sense, antisense = get_dual_strand("AUGCAUGC")
    assert len(sense) == len(antisense)
    assert "U" not in sense  # should be DNA
    assert sense == "ATGCATGC"
    assert antisense == reverse_complement_dna(sense)


def test_find_pam_sites_NGG():
    # Place a known NGG PAM preceded by a 20-nt spacer
    #   spacer (20 nt)            PAM
    seq = "ACGTACGTACGTACGTACGT" + "AGG" + "TTTTTT"
    matches = find_pam_sites(seq, "NGG", "3prime", 20)
    assert len(matches) >= 1
    sense_hits = [m for m in matches if m.strand == "sense" and m.pam_seq == "AGG"]
    assert len(sense_hits) >= 1
    hit = sense_hits[0]
    assert hit.spacer_seq == "ACGTACGTACGTACGTACGT"


def test_find_pam_sites_no_match():
    # Sequence with no NGG anywhere (all As)
    seq = "A" * 40
    matches = find_pam_sites(seq, "NGG", "3prime", 20)
    # Sense strand has no G at all, antisense is all T -> no NGG
    assert matches == []


def test_gc_content_all_gc():
    assert gc_content("GCGC") == 1.0


def test_gc_content_no_gc():
    assert gc_content("AATT") == 0.0


def test_gc_content_empty():
    assert gc_content("") == 0.0


def test_has_homopolymer_true():
    assert has_homopolymer("AAAAAATG", 4) is True


def test_has_homopolymer_false():
    assert has_homopolymer("ATGATG", 4) is False


def test_count_mismatches_perfect():
    score = count_mismatches_weighted("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT")
    assert score == 0.0


def test_count_mismatches_single():
    # Mismatch at position 0 (seed region) should use weight 2.0
    guide  = "ACGTACGTACGTACGTACGT"
    target = "XCGTACGTACGTACGTACGT"
    score = count_mismatches_weighted(guide, target)
    assert score == 2.0  # position 0 is in seed (weight 2.0)
