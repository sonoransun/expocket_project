"""Tests for viewer.genome_editing.guide_designer."""

from __future__ import annotations

from viewer.genome_editing.guide_designer import (
    design_cas9_guide,
    design_cas12a_guide,
    score_guide,
)
from viewer.genome_editing.tools import TargetSite


def _make_site(tool: str, spacer: str, pam: str = "AGG", strand: str = "sense") -> TargetSite:
    """Helper to build a minimal TargetSite for testing."""
    return TargetSite(
        tool=tool,
        variant_id="A_AAA",
        strand=strand,
        spacer_start=0,
        spacer_end=len(spacer),
        spacer_seq=spacer,
        pam_seq=pam,
        cut_positions=(17, 17),
        gc_content=0.5,
        seed_score=0.8,
        binding_score=0.8,
    )


def test_design_cas9_guide():
    site = _make_site("cas9", "ACGTACGTACGTACGTACGT")
    spacer, scaffold = design_cas9_guide(site)
    assert spacer == "ACGTACGTACGTACGTACGT"
    assert len(scaffold) > 0
    # scaffold is the canonical tracrRNA scaffold
    assert "GUUUUAGAGC" in scaffold or "GTTTTAGAGC" in scaffold


def test_design_cas12a_guide():
    site = _make_site("cas12a", "ACGTACGTACGTACGTACGTG", pam="TTTC")
    direct_repeat, spacer = design_cas12a_guide(site)
    assert spacer == "ACGTACGTACGTACGTACGTG"
    assert len(direct_repeat) > 0
    assert direct_repeat == "AATTTCTACTCTTGTAGAT"


def test_score_guide_ideal():
    # 50% GC, no homopolymer, moderate seed GC -> high score
    spacer = "ACGTACGTACGTACGTACGT"  # 50% GC
    s = score_guide(spacer)
    assert s > 0.7


def test_score_guide_extreme_gc():
    # All G -> GC=100%, also has GGGG homopolymer -> low score
    spacer = "G" * 20
    s = score_guide(spacer)
    assert s < 0.5
