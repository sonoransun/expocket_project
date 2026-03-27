"""Generate guide sequences and binding domain descriptions for each editing tool."""

from __future__ import annotations

from viewer.genome_editing.tools import TargetSite
from viewer.genome_editing.sequence_utils import (
    reverse_complement_dna,
    gc_content,
    has_homopolymer,
)

# ---------------------------------------------------------------------------
# Scaffold / repeat sequences
# ---------------------------------------------------------------------------
_CAS9_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
)
_CAS12A_DIRECT_REPEAT = "AATTTCTACTCTTGTAGAT"

# RVD table (for TALEN arm string generation)
_RVD: dict[str, str] = {"A": "NI", "C": "HD", "G": "NN", "T": "NG"}


def _guide_for_strand(site: TargetSite) -> str:
    """Return the spacer written 5'→3' for use in the expression vector."""
    if site.strand == "sense":
        return site.spacer_seq
    # Antisense site: the expressed guide matches the antisense strand spacer
    return site.spacer_seq


def design_cas9_guide(site: TargetSite) -> tuple[str, str]:
    """Return (spacer_20nt_DNA, sgRNA_scaffold)."""
    return _guide_for_strand(site), _CAS9_SCAFFOLD


def design_cas12a_guide(site: TargetSite) -> tuple[str, str]:
    """Return (direct_repeat_DNA, spacer_21nt_DNA) for Cas12a crRNA."""
    return _CAS12A_DIRECT_REPEAT, _guide_for_strand(site)


def design_casclover_pair(site: TargetSite) -> tuple[str, str, str, str]:
    """Return (guide1_spacer, guide1_scaffold, guide2_spacer, guide2_scaffold)."""
    g1 = _guide_for_strand(site)
    s1 = _CAS9_SCAFFOLD
    if site.paired_site is not None:
        g2 = _guide_for_strand(site.paired_site)
    else:
        g2 = ""
    s2 = _CAS9_SCAFFOLD
    return g1, s1, g2, s2


def design_nicer_pair(site: TargetSite) -> tuple[str, str, str, str]:
    """Return (guide1_spacer, guide1_scaffold, guide2_spacer, guide2_scaffold) for NICER."""
    return design_casclover_pair(site)


def design_talen_arms(site: TargetSite) -> dict:
    """Return a dict with left and right TALEN arm RVD strings and sequences."""
    arms = site.binding_arms
    if not arms:
        return {}
    left_rvd = arms.get("left_rvd", "")
    right_rvd = arms.get("right_rvd", "")
    return {
        "left_rvd_string": left_rvd,
        "right_rvd_string": right_rvd,
        "left_dna": arms.get("left_seq", ""),
        "right_dna": arms.get("right_seq", ""),
        "spacer": arms.get("spacer_seq", ""),
        "spacer_len": arms.get("spacer_len", 0),
        "left_arm_len": arms.get("left_len", 0),
        "right_arm_len": arms.get("right_len", 0),
        "note": (
            f"Left arm ({arms.get('left_len', 0)} aa RVD repeats): {left_rvd}\n"
            f"Right arm ({arms.get('right_len', 0)} aa RVD repeats): {right_rvd}\n"
            f"Spacer ({arms.get('spacer_len', 0)} bp) — FokI cuts in spacer"
        ),
    }


def design_zfn_fingers(site: TargetSite) -> dict:
    """Return a dict describing left and right zinc finger arrays."""
    arms = site.binding_arms
    if not arms:
        return {}
    return {
        "left_fingers": arms.get("left_fingers", []),
        "right_fingers": arms.get("right_fingers", []),
        "left_dna": arms.get("left_seq", ""),
        "right_dna": arms.get("right_seq", ""),
        "spacer_len": arms.get("spacer_len", 0),
        "spacer_seq": arms.get("spacer_seq", ""),
        "note": (
            f"Left ZF array: "
            f"{' | '.join(f['triplet'] for f in arms.get('left_fingers', []))}\n"
            f"Right ZF array: "
            f"{' | '.join(f['triplet'] for f in arms.get('right_fingers', []))}"
        ),
    }


def design_retron_template(site: TargetSite, edit_seq: str, flank: int = 20) -> str:
    """Build an EC86-style retron msd template for the edit at site.

    Format: 5'-NCGCAATG + left_flank(flank nt) + edit_seq + right_flank(flank nt) + CATTGCGN-3'
    The flanking sequences come from the sense-strand DNA context around the edit position.
    """
    from viewer.genome_editing.sequence_utils import get_dual_strand
    # We need the full sequence; use spacer_seq context heuristic
    # The site already carries spacer_start (= edit position for retron)
    left_context = f"{'N' * max(0, flank - site.spacer_start)}"
    right_context = "N" * flank  # will be filled by edit_engine when full seq is available
    return f"NCGCAATG{left_context}{edit_seq}{right_context}CATTGCGN"


def score_guide(spacer: str) -> float:
    """0-1 quality score: GC content, no TTTT, homopolymer penalty, seed GC balance."""
    gc = gc_content(spacer)
    if 0.4 <= gc <= 0.7:
        gc_score = 1.0
    elif gc < 0.4:
        gc_score = gc / 0.4
    else:
        gc_score = max(0.0, 1.0 - (gc - 0.7) / 0.3)

    penalty = 0.0
    if "TTTT" in spacer.upper():
        penalty += 0.4
    elif has_homopolymer(spacer, 4):
        penalty += 0.2

    seed = spacer[-8:] if len(spacer) >= 8 else spacer
    seed_gc = gc_content(seed)
    if 0.35 <= seed_gc <= 0.65:
        seed_score = 1.0
    else:
        seed_score = 0.5

    return max(0.0, min(1.0, (gc_score * 0.5 + seed_score * 0.5) - penalty))


def guide_summary(site: TargetSite) -> str:
    """Return a short human-readable guide or binding domain description."""
    tool = site.tool
    if tool == "cas9":
        spacer, scaffold = design_cas9_guide(site)
        return f"sgRNA spacer (5'→3'): {spacer}"
    if tool == "cas12a":
        dr, spacer = design_cas12a_guide(site)
        return f"crRNA spacer (5'→3'): {spacer}  | direct-repeat: {dr[:12]}…"
    if tool in ("casclover", "nicer"):
        g1, _, g2, _ = design_casclover_pair(site)
        return f"Guide1: {g1}  | Guide2: {g2 or 'N/A'}"
    if tool == "talen":
        d = design_talen_arms(site)
        la = d.get("left_dna", "")
        ra = d.get("right_dna", "")
        return f"Left arm: {la[:18]}… | Right arm: {ra[:18]}…"
    if tool == "zfn":
        d = design_zfn_fingers(site)
        lt = [f['triplet'] for f in d.get("left_fingers", [])]
        rt = [f['triplet'] for f in d.get("right_fingers", [])]
        return f"Left: {'-'.join(lt)} | Right: {'-'.join(rt)}"
    if tool == "retron":
        return f"Edit position: {site.spacer_start} (no guide required — template-directed)"
    return ""
