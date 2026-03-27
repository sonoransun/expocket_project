"""RNA/DNA conversion, reverse complement, and PAM-site scanning utilities."""

from __future__ import annotations

import re
from dataclasses import dataclass

_RNA_TO_DNA = str.maketrans("UuAaGgCc", "TtAaGgCc")
_DNA_TO_RNA = str.maketrans("TtAaGgCc", "UuAaGgCc")
_DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

IUPAC_REGEX: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "N": "[ACGT]", "R": "[AG]", "Y": "[CT]",
    "S": "[GC]", "W": "[AT]", "K": "[GT]",
    "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]",
}


def rna_to_dna(seq: str) -> str:
    """Replace U with T (case-preserving)."""
    return seq.translate(_RNA_TO_DNA)


def dna_to_rna(seq: str) -> str:
    """Replace T with U (case-preserving)."""
    return seq.translate(_DNA_TO_RNA)


def reverse_complement_dna(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.upper().translate(_DNA_COMP)[::-1]


def reverse_complement_rna(seq: str) -> str:
    """Return the reverse complement of an RNA sequence."""
    _rna_comp = str.maketrans("ACGUacgu", "UGCAugca")
    return seq.upper().translate(_rna_comp)[::-1]


def get_dual_strand(rna_seq: str) -> tuple[str, str]:
    """Convert RNA hairpin to sense and antisense DNA strands.

    Returns (sense_dna, antisense_dna) in 5'→3' orientation.
    sense_dna:     direct U→T of rna_seq (5'→3')
    antisense_dna: reverse complement (5'→3')
    """
    sense = rna_to_dna(rna_seq).upper()
    antisense = reverse_complement_dna(sense)
    return sense, antisense


def _pam_to_regex(pam: str) -> str:
    """Convert an IUPAC PAM string to a compiled regex pattern."""
    return "".join(IUPAC_REGEX.get(c.upper(), c) for c in pam)


@dataclass
class PamMatch:
    """A single PAM + spacer match on a given strand."""
    strand: str           # "sense" or "antisense"
    pam_start: int        # in sense-strand coordinates (0-based)
    pam_end: int
    spacer_start: int     # in sense-strand coordinates (0-based)
    spacer_end: int
    spacer_seq: str       # spacer DNA, always written 5'→3'
    pam_seq: str          # actual matched PAM bases


def find_pam_sites(
    sense_dna: str,
    pam_pattern: str,
    pam_position: str,
    spacer_len: int,
) -> list[PamMatch]:
    """Find all PAM sites on both strands of a DNA sequence.

    Args:
        sense_dna:    sense (non-template) DNA strand, 5'→3'
        pam_pattern:  IUPAC PAM string (e.g. "NGG", "TTTN")
        pam_position: "3prime" — PAM is 3' of spacer (Cas9 style)
                      "5prime" — PAM is 5' of spacer (Cas12a style)
        spacer_len:   guide spacer length in nt

    Returns list of PamMatch (all positions in sense-strand coordinates).
    """
    pam_re = re.compile(_pam_to_regex(pam_pattern), re.IGNORECASE)
    pam_len = len(pam_pattern)
    seq_len = len(sense_dna)
    results: list[PamMatch] = []

    antisense = reverse_complement_dna(sense_dna)

    for strand_label, seq in (("sense", sense_dna), ("antisense", antisense)):
        for m in pam_re.finditer(seq):
            p_start = m.start()
            p_end = m.end()

            if pam_position == "3prime":
                # spacer is immediately 5' of PAM
                sp_end = p_start
                sp_start = p_start - spacer_len
                if sp_start < 0:
                    continue
            else:
                # 5prime: spacer is immediately 3' of PAM
                sp_start = p_end
                sp_end = p_end + spacer_len
                if sp_end > len(seq):
                    continue

            spacer_seq = seq[sp_start:sp_end]

            # Convert positions to sense-strand coordinates
            if strand_label == "sense":
                sense_p_start = p_start
                sense_p_end = p_end
                sense_sp_start = sp_start
                sense_sp_end = sp_end
            else:
                # Antisense positions need to be mirrored to sense coords
                # antisense index i → sense index (seq_len - 1 - i)
                sense_p_start = seq_len - p_end
                sense_p_end = seq_len - p_start
                sense_sp_start = seq_len - sp_end
                sense_sp_end = seq_len - sp_start

            results.append(PamMatch(
                strand=strand_label,
                pam_start=sense_p_start,
                pam_end=sense_p_end,
                spacer_start=sense_sp_start,
                spacer_end=sense_sp_end,
                spacer_seq=spacer_seq.upper(),
                pam_seq=m.group(0).upper(),
            ))

    return results


def gc_content(seq: str) -> float:
    """Return GC fraction (0-1) of a nucleotide sequence."""
    seq_u = seq.upper()
    gc = seq_u.count("G") + seq_u.count("C")
    return gc / len(seq) if seq else 0.0


def has_homopolymer(seq: str, length: int = 4) -> bool:
    """Return True if seq contains a homopolymer run >= length."""
    for nt in "ACGT":
        if nt * length in seq.upper():
            return True
    return False


def count_mismatches_weighted(
    guide: str,
    target: str,
    weights: list[float] | None = None,
) -> float:
    """Compute weighted mismatch score between guide and target (same length).

    Position 0 = seed-proximal (PAM-proximal) end.
    Defaults: positions 0-7 (seed) weight 2.0, rest weight 0.5.
    Returns total weighted mismatches.
    """
    n = min(len(guide), len(target))
    if weights is None:
        weights = [2.0 if i < 8 else 0.5 for i in range(n)]
    score = 0.0
    for i in range(n):
        if guide[i].upper() != target[i].upper():
            score += weights[i] if i < len(weights) else 0.5
    return score
