"""Find valid target sites for each gene editing tool by scanning both DNA strands."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from viewer.genome_editing.sequence_utils import (
    gc_content,
    get_dual_strand,
    has_homopolymer,
    find_pam_sites,
    reverse_complement_dna,
    count_mismatches_weighted,
)
from viewer.genome_editing.tools import TOOLS, TargetSite

if TYPE_CHECKING:
    from viewer.data.schema import EnrichedVariantDataset

# ---------------------------------------------------------------------------
# TALEN RVD table
# ---------------------------------------------------------------------------
_RVD: dict[str, str] = {
    "A": "NI", "C": "HD", "G": "NN", "T": "NG",
}

# ---------------------------------------------------------------------------
# ZFN simplified triplet recognition table (Zif268-derived + common designs)
# ---------------------------------------------------------------------------
_ZFN_TRIPLETS: dict[str, str] = {
    "GCG": "C2H2-Zif268F2", "GCC": "C2H2-Zif268F1", "GCT": "C2H2-type-I",
    "GCA": "C2H2-type-II",  "GAG": "C2H2-EGR1f3",  "GAA": "C2H2-EGR1f2",
    "GGG": "C2H2-WT-1",     "GGC": "C2H2-GGC",     "GAC": "C2H2-GAC",
    "GTG": "C2H2-GTG",      "GTC": "C2H2-GTC",      "GTA": "C2H2-GTA",
    "ACG": "C2H2-ACG",      "ACC": "C2H2-ACC",      "ACA": "C2H2-ACA",
    "TCG": "C2H2-TCG",      "TCC": "C2H2-TCC",      "TCA": "C2H2-TCA",
    "CAG": "C2H2-CAG",      "CGG": "C2H2-CGG",
}

# EC86 retron RT recognition handles (simplified)
_EC86_LEFT_HANDLE = "NCGCAATG"
_EC86_RIGHT_HANDLE = "CATTGCGN"


def _score_guide(spacer: str) -> float:
    """Composite 0-1 quality score for a 20-21 nt DNA guide spacer."""
    gc = gc_content(spacer)
    # GC content: ideal 40-70%
    if 0.4 <= gc <= 0.7:
        gc_score = 1.0
    elif gc < 0.4:
        gc_score = gc / 0.4
    else:
        gc_score = max(0.0, 1.0 - (gc - 0.7) / 0.3)

    # Penalise TTTT (Pol III terminator) and homopolymer > 4
    penalty = 0.0
    if "TTTT" in spacer.upper():
        penalty += 0.4
    elif has_homopolymer(spacer, 4):
        penalty += 0.2

    # Seed GC balance (positions 1-8 from PAM side = last 8 of spacer)
    seed = spacer[-8:] if len(spacer) >= 8 else spacer
    seed_gc = gc_content(seed)
    if 0.35 <= seed_gc <= 0.65:
        seed_score = 1.0
    else:
        seed_score = 0.5

    return max(0.0, min(1.0, (gc_score * 0.5 + seed_score * 0.5) - penalty))


def _find_cas9_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    matches = find_pam_sites(sense_dna, "NGG", "3prime", 20)
    sites = []
    for m in matches:
        gc = gc_content(m.spacer_seq)
        score = _score_guide(m.spacer_seq)
        # Cut between pos 17-18 of spacer (3 bp from PAM)
        if m.strand == "sense":
            cut_sense = m.spacer_start + 17
            cut_anti = cut_sense       # blunt
        else:
            cut_sense = m.spacer_end - 17
            cut_anti = cut_sense
        sites.append(TargetSite(
            tool="cas9",
            variant_id=variant_id,
            strand=m.strand,
            spacer_start=m.spacer_start,
            spacer_end=m.spacer_end,
            spacer_seq=m.spacer_seq,
            pam_seq=m.pam_seq,
            cut_positions=(cut_sense, cut_anti),
            gc_content=gc,
            seed_score=score,
            binding_score=score,
        ))
    return sorted(sites, key=lambda s: -s.binding_score)


def _find_cas12a_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    matches = find_pam_sites(sense_dna, "TTTN", "5prime", 21)
    sites = []
    for m in matches:
        gc = gc_content(m.spacer_seq)
        score = _score_guide(m.spacer_seq)
        # Staggered cut: 18/23 nt downstream of PAM → 5-nt 5' overhang
        if m.strand == "sense":
            cut_sense = m.spacer_start + 18
            cut_anti = m.spacer_start + 23
        else:
            cut_sense = m.spacer_end - 18
            cut_anti = m.spacer_end - 23
        sites.append(TargetSite(
            tool="cas12a",
            variant_id=variant_id,
            strand=m.strand,
            spacer_start=m.spacer_start,
            spacer_end=m.spacer_end,
            spacer_seq=m.spacer_seq,
            pam_seq=m.pam_seq,
            cut_positions=(cut_sense, cut_anti),
            gc_content=gc,
            seed_score=score,
            binding_score=score,
        ))
    return sorted(sites, key=lambda s: -s.binding_score)


def _find_paired_cas9n_sites(
    sense_dna: str,
    variant_id: str,
    tool_name: str,
    window_min: int,
    window_max: int,
) -> list[TargetSite]:
    """Find paired Cas9n (D10A) sites for Cas-CLOVER or NICER.

    Pairs: one guide on sense strand + one guide on antisense strand,
    both NGG PAMs pointing outward (PAM-out orientation),
    with spacer_end gap inside window_min..window_max.
    """
    sense_matches = [m for m in find_pam_sites(sense_dna, "NGG", "3prime", 20)
                     if m.strand == "sense"]
    anti_matches = [m for m in find_pam_sites(sense_dna, "NGG", "3prime", 20)
                    if m.strand == "antisense"]

    sites: list[TargetSite] = []
    for s in sense_matches:
        for a in anti_matches:
            # PAM-out: sense guide PAM is 3' of spacer (to the right),
            # antisense guide PAM is also 3' of its spacer on the antisense strand
            # → on sense coords, antisense spacer_start < sense spacer_end
            # Gap = a.spacer_start - s.spacer_end  (sense coords)
            gap = a.spacer_start - s.spacer_end
            if window_min <= gap <= window_max:
                score_s = _score_guide(s.spacer_seq)
                score_a = _score_guide(a.spacer_seq)
                composite = (score_s + score_a) / 2

                # Nick positions (D10A nicks RuvC strand = sense strand nicked by sense guide,
                # antisense strand nicked by antisense guide)
                nick1 = s.spacer_start + 17  # sense nick
                nick2 = a.spacer_start + 17  # antisense nick (in sense coords)

                sense_site = TargetSite(
                    tool=tool_name,
                    variant_id=variant_id,
                    strand="sense",
                    spacer_start=s.spacer_start,
                    spacer_end=s.spacer_end,
                    spacer_seq=s.spacer_seq,
                    pam_seq=s.pam_seq,
                    cut_positions=(nick1, nick2),
                    gc_content=gc_content(s.spacer_seq),
                    seed_score=score_s,
                    binding_score=composite,
                )
                anti_site = TargetSite(
                    tool=tool_name,
                    variant_id=variant_id,
                    strand="antisense",
                    spacer_start=a.spacer_start,
                    spacer_end=a.spacer_end,
                    spacer_seq=a.spacer_seq,
                    pam_seq=a.pam_seq,
                    cut_positions=(nick1, nick2),
                    gc_content=gc_content(a.spacer_seq),
                    seed_score=score_a,
                    binding_score=composite,
                )
                sense_site.paired_site = anti_site
                sites.append(sense_site)

    return sorted(sites, key=lambda s: -s.binding_score)[:20]


def _find_talen_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    """Find TALEN binding sites: left_arm + spacer + right_arm."""
    seq = sense_dna.upper()
    seq_len = len(seq)
    sites: list[TargetSite] = []

    left_lens = range(14, 21)
    spacer_lens = range(12, 21)
    right_lens = range(14, 21)

    for sp_len in spacer_lens:
        for la_len in left_lens:
            for ra_len in right_lens:
                total = la_len + sp_len + ra_len
                for start in range(0, seq_len - total + 1):
                    left_seq = seq[start: start + la_len]
                    spacer_seq = seq[start + la_len: start + la_len + sp_len]
                    right_seq = seq[start + la_len + sp_len: start + total]

                    # Reject T-rich stretches (TALE binding dislikes homopolymers)
                    if has_homopolymer(left_seq, 5) or has_homopolymer(right_seq, 5):
                        continue

                    # RVD codes
                    left_rvd = "-".join(_RVD.get(nt, "NN") for nt in left_seq)
                    right_rvd = "-".join(_RVD.get(nt, "NN") for nt in right_seq)

                    # Score: prefer moderate GC, no homopolymers
                    gc_l = gc_content(left_seq)
                    gc_r = gc_content(right_seq)
                    score = 1.0 - 0.5 * (abs(gc_l - 0.5) + abs(gc_r - 0.5))

                    # FokI cuts in spacer midpoint
                    cut = start + la_len + sp_len // 2
                    sites.append(TargetSite(
                        tool="talen",
                        variant_id=variant_id,
                        strand="sense",
                        spacer_start=start,
                        spacer_end=start + total,
                        spacer_seq=spacer_seq,
                        pam_seq="",
                        cut_positions=(cut, cut + 1),
                        gc_content=(gc_l + gc_r) / 2,
                        seed_score=score,
                        binding_score=score,
                        binding_arms={
                            "left_seq": left_seq,
                            "left_rvd": left_rvd,
                            "left_len": la_len,
                            "right_seq": right_seq,
                            "right_rvd": right_rvd,
                            "right_len": ra_len,
                            "spacer_len": sp_len,
                            "spacer_seq": spacer_seq,
                        },
                    ))

    # Deduplicate by spacer midpoint, keep best score
    best: dict[int, TargetSite] = {}
    for s in sites:
        mid = (s.spacer_start + s.spacer_end) // 2
        if mid not in best or s.binding_score > best[mid].binding_score:
            best[mid] = s
    return sorted(best.values(), key=lambda s: -s.binding_score)[:20]


def _find_zfn_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    """Find ZFN binding sites: left arm (3-6 triplets) + spacer + right arm."""
    seq = sense_dna.upper()
    seq_len = len(seq)
    sites: list[TargetSite] = []

    finger_counts = range(3, 7)   # 3-6 fingers per arm
    spacer_lens = range(5, 8)      # 5-7 bp spacer

    for n_left in finger_counts:
        for n_right in finger_counts:
            la_len = n_left * 3
            ra_len = n_right * 3
            for sp_len in spacer_lens:
                total = la_len + sp_len + ra_len
                for start in range(0, seq_len - total + 1):
                    left_seq = seq[start: start + la_len]
                    spacer_seq = seq[start + la_len: start + la_len + sp_len]
                    right_seq = seq[start + la_len + sp_len: start + total]

                    # Parse left arm triplets (read 3'→5' for ZF convention)
                    left_triplets = [left_seq[i:i+3] for i in range(0, la_len, 3)]
                    right_triplets = [right_seq[i:i+3] for i in range(0, ra_len, 3)]

                    # Count recognizable triplets
                    rec_left = sum(1 for t in left_triplets if t in _ZFN_TRIPLETS)
                    rec_right = sum(1 for t in right_triplets if t in _ZFN_TRIPLETS)
                    if rec_left < n_left - 1 or rec_right < n_right - 1:
                        continue

                    score = (rec_left + rec_right) / (n_left + n_right)
                    cut = start + la_len + sp_len // 2

                    left_finger_info = [
                        {"triplet": t, "zf_type": _ZFN_TRIPLETS.get(t, "unknown"),
                         "position": start + i * 3}
                        for i, t in enumerate(left_triplets)
                    ]
                    right_finger_info = [
                        {"triplet": t, "zf_type": _ZFN_TRIPLETS.get(t, "unknown"),
                         "position": start + la_len + sp_len + i * 3}
                        for i, t in enumerate(right_triplets)
                    ]

                    sites.append(TargetSite(
                        tool="zfn",
                        variant_id=variant_id,
                        strand="sense",
                        spacer_start=start,
                        spacer_end=start + total,
                        spacer_seq=spacer_seq,
                        pam_seq="",
                        cut_positions=(cut, cut - 4),  # FokI 4-nt 3' overhang
                        gc_content=gc_content(left_seq + right_seq),
                        seed_score=score,
                        binding_score=score,
                        binding_arms={
                            "left_seq": left_seq,
                            "left_fingers": left_finger_info,
                            "right_seq": right_seq,
                            "right_fingers": right_finger_info,
                            "spacer_len": sp_len,
                            "spacer_seq": spacer_seq,
                        },
                    ))

    best: dict[int, TargetSite] = {}
    for s in sites:
        mid = (s.spacer_start + s.spacer_end) // 2
        if mid not in best or s.binding_score > best[mid].binding_score:
            best[mid] = s
    return sorted(best.values(), key=lambda s: -s.binding_score)[:20]


def _find_retron_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    """EC86 retron: one site per position (no PAM constraint)."""
    seq = sense_dna.upper()
    sites: list[TargetSite] = []
    flank = 20  # nt of flanking sequence for template
    for pos in range(flank, len(seq) - flank):
        left_flank = seq[pos - flank: pos]
        right_flank = seq[pos: pos + flank]
        # Score based on flank uniqueness (simple: reward non-repetitive flanks)
        score = 0.8 - 0.3 * (has_homopolymer(left_flank + right_flank, 5))
        sites.append(TargetSite(
            tool="retron",
            variant_id=variant_id,
            strand="sense",
            spacer_start=pos,
            spacer_end=pos + 1,
            spacer_seq=seq[pos],
            pam_seq="",
            cut_positions=(pos, pos),
            gc_content=gc_content(left_flank + right_flank),
            seed_score=score,
            binding_score=score,
        ))
    return sites  # all positions returned for this tool


def _find_nicer_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    return _find_paired_cas9n_sites(sense_dna, variant_id, "nicer", 5, 25)


def _find_casclover_sites(sense_dna: str, variant_id: str) -> list[TargetSite]:
    return _find_paired_cas9n_sites(sense_dna, variant_id, "casclover", 25, 50)


_TOOL_FINDERS = {
    "cas9": _find_cas9_sites,
    "cas12a": _find_cas12a_sites,
    "casclover": _find_casclover_sites,
    "talen": _find_talen_sites,
    "zfn": _find_zfn_sites,
    "retron": _find_retron_sites,
    "nicer": _find_nicer_sites,
}


class TargetFinder:
    """Facade: finds valid editing sites for any tool on any variant."""

    def find_sites(
        self,
        variant_id: str,
        tool_name: str,
        dataset: "EnrichedVariantDataset",
    ) -> list[TargetSite]:
        """Return ranked target sites for tool_name on variant_id."""
        variant = dataset.get_variant(variant_id)
        if variant is None:
            return []
        sense_dna, _ = get_dual_strand(variant.pre_mirna_sequence)
        finder = _TOOL_FINDERS.get(tool_name)
        if finder is None:
            return []
        return finder(sense_dna, variant_id)

    def find_all_tools(
        self,
        variant_id: str,
        dataset: "EnrichedVariantDataset",
        tool_names: list[str] | None = None,
    ) -> dict[str, list[TargetSite]]:
        """Return sites for all (or a subset of) tools."""
        names = tool_names or list(_TOOL_FINDERS.keys())
        return {name: self.find_sites(variant_id, name, dataset) for name in names}
