"""Sequence editing simulation: SNP, insertion, deletion, HDR template generation."""

from __future__ import annotations

from viewer.genome_editing.sequence_utils import (
    rna_to_dna,
    dna_to_rna,
    reverse_complement_dna,
    get_dual_strand,
)

# PAM-silencing mutations (to prevent re-cutting after HDR)
_PAM_SILENCE: dict[str, str] = {
    "NGG": {"GG": "GA"},   # mutate 2nd G of NGG to A
    "TTTN": {"TTT": "TCT"},  # mutate middle T of TTTN
}


def apply_snp(dna: str, position: int, new_base: str) -> str:
    """Replace one base at position with new_base (DNA)."""
    if 0 <= position < len(dna):
        return dna[:position] + new_base.upper() + dna[position + 1:]
    return dna


def apply_insertion(dna: str, position: int, insert: str) -> str:
    """Insert bases after position."""
    return dna[:position] + insert.upper() + dna[position:]


def apply_deletion(dna: str, position: int, length: int) -> str:
    """Delete `length` bases starting at position."""
    end = min(position + length, len(dna))
    return dna[:position] + dna[end:]


def _silence_pam(dna: str, cut_pos: int, pam_pattern: str | None) -> str:
    """Silently mutate PAM near cut_pos to prevent Cas re-cutting after HDR."""
    if pam_pattern is None or cut_pos < 0 or cut_pos >= len(dna):
        return dna
    # Cas9 NGG: PAM is 3 nt after cut_pos (approx)
    window_start = max(0, cut_pos)
    window_end = min(len(dna), cut_pos + 6)
    window = dna[window_start:window_end]
    for pam_core, replacement in _PAM_SILENCE.get(pam_pattern, {}).items():
        if pam_core in window:
            new_window = window.replace(pam_core, replacement, 1)
            return dna[:window_start] + new_window + dna[window_end:]
    return dna


def generate_hdr_template(
    sense_dna: str,
    cut_pos: int,
    edit_position: int,
    new_base: str,
    flank: int = 50,
    pam_pattern: str | None = None,
) -> str:
    """Generate a ~100 nt ssODN HDR template.

    Applies the edit at edit_position, silences the PAM near cut_pos,
    and returns flank nt left + edit nt + flank nt right.
    """
    seq = sense_dna.upper()
    # Apply edit
    edited = apply_snp(seq, edit_position, new_base)
    # Silence PAM to prevent re-cutting
    edited = _silence_pam(edited, cut_pos, pam_pattern)
    # Trim to template window around cut
    left_start = max(0, cut_pos - flank)
    right_end = min(len(edited), cut_pos + flank)
    return edited[left_start:right_end]


def generate_insertion_hdr_template(
    sense_dna: str,
    cut_pos: int,
    insert_position: int,
    insert_seq: str,
    flank: int = 50,
    pam_pattern: str | None = None,
) -> str:
    """HDR template for an insertion."""
    seq = sense_dna.upper()
    edited = apply_insertion(seq, insert_position, insert_seq)
    edited = _silence_pam(edited, cut_pos, pam_pattern)
    left_start = max(0, cut_pos - flank)
    right_end = min(len(edited), cut_pos + flank + len(insert_seq))
    return edited[left_start:right_end]


def generate_deletion_hdr_template(
    sense_dna: str,
    cut_pos: int,
    del_position: int,
    del_length: int,
    flank: int = 50,
    pam_pattern: str | None = None,
) -> str:
    """HDR template for a deletion."""
    seq = sense_dna.upper()
    edited = apply_deletion(seq, del_position, del_length)
    edited = _silence_pam(edited, cut_pos, pam_pattern)
    left_start = max(0, cut_pos - flank)
    right_end = min(len(edited), cut_pos + flank)
    return edited[left_start:right_end]


def predict_nhej_outcomes(sense_dna: str, cut_pos: int, n: int = 5) -> list[str]:
    """Return n likely NHEJ indel outcomes (simplified).

    Based on empirical NHEJ rules:
    - 1-bp deletion at cut is most common
    - 1-bp insertion (duplication of adjacent base) second
    - Larger deletions less frequent
    """
    seq = sense_dna.upper()
    outcomes: list[str] = []
    if n >= 1:
        outcomes.append(apply_deletion(seq, cut_pos, 1))      # 1-bp del
    if n >= 2:
        dup = seq[cut_pos] if cut_pos < len(seq) else "N"
        outcomes.append(apply_insertion(seq, cut_pos, dup))   # 1-bp ins (dup)
    if n >= 3:
        outcomes.append(apply_deletion(seq, max(0, cut_pos - 1), 2))  # 2-bp del
    if n >= 4:
        outcomes.append(apply_deletion(seq, cut_pos, 3))      # 3-bp del
    if n >= 5:
        outcomes.append(apply_deletion(seq, max(0, cut_pos - 2), 4))  # 4-bp del
    return outcomes[:n]


def design_retron_msd(
    sense_dna: str,
    edit_position: int,
    edit_base: str,
    flank: int = 20,
) -> str:
    """Build an EC86 retron msd template.

    Format: 5'-NCGCAATG + left_flank(flank) + edit + right_flank(flank) + CATTGCGN-3'
    """
    seq = sense_dna.upper()
    left_start = max(0, edit_position - flank)
    right_end = min(len(seq), edit_position + 1 + flank)
    left_flank = seq[left_start:edit_position]
    right_flank = seq[edit_position + 1:right_end]
    return f"NCGCAATG{left_flank}{edit_base.upper()}{right_flank}CATTGCGN"


def rna_sequence_after_snp(original_rna: str, position: int, new_base_dna: str) -> str:
    """Apply a SNP edit (in DNA space) and return the resulting RNA sequence."""
    sense_dna = rna_to_dna(original_rna).upper()
    edited_dna = apply_snp(sense_dna, position, new_base_dna)
    return dna_to_rna(edited_dna)


def rna_sequence_after_insertion(original_rna: str, position: int, insert_dna: str) -> str:
    sense_dna = rna_to_dna(original_rna).upper()
    edited_dna = apply_insertion(sense_dna, position, insert_dna)
    return dna_to_rna(edited_dna)


def rna_sequence_after_deletion(original_rna: str, position: int, length: int) -> str:
    sense_dna = rna_to_dna(original_rna).upper()
    edited_dna = apply_deletion(sense_dna, position, length)
    return dna_to_rna(edited_dna)


class EditEngine:
    """Stateless façade for sequence editing operations."""

    def apply_edit(
        self,
        original_rna: str,
        edit_type: str,
        edit_position: int,
        edit_sequence: str,
        cut_pos: int = -1,
        pam_pattern: str | None = None,
    ) -> tuple[str, str]:
        """Apply an edit to original_rna.

        Returns (modified_rna, hdr_template).
        """
        sense_dna = rna_to_dna(original_rna).upper()

        if edit_type == "snp":
            modified_rna = rna_sequence_after_snp(original_rna, edit_position, edit_sequence)
            hdr = generate_hdr_template(
                sense_dna, cut_pos if cut_pos >= 0 else edit_position,
                edit_position, edit_sequence, pam_pattern=pam_pattern,
            )
        elif edit_type == "insertion":
            modified_rna = rna_sequence_after_insertion(original_rna, edit_position, edit_sequence)
            hdr = generate_insertion_hdr_template(
                sense_dna, cut_pos if cut_pos >= 0 else edit_position,
                edit_position, edit_sequence, pam_pattern=pam_pattern,
            )
        elif edit_type == "deletion":
            length = len(edit_sequence) if edit_sequence else 1
            modified_rna = rna_sequence_after_deletion(original_rna, edit_position, length)
            hdr = generate_deletion_hdr_template(
                sense_dna, cut_pos if cut_pos >= 0 else edit_position,
                edit_position, length, pam_pattern=pam_pattern,
            )
        else:
            modified_rna = original_rna
            hdr = ""

        return modified_rna, hdr

    def retron_template(
        self,
        original_rna: str,
        edit_position: int,
        edit_base_dna: str,
        flank: int = 20,
    ) -> str:
        sense_dna = rna_to_dna(original_rna).upper()
        return design_retron_msd(sense_dna, edit_position, edit_base_dna, flank)

    def nhej_outcomes(self, original_rna: str, cut_pos: int, n: int = 5) -> list[str]:
        sense_dna = rna_to_dna(original_rna).upper()
        dna_outcomes = predict_nhej_outcomes(sense_dna, cut_pos, n)
        return [dna_to_rna(d) for d in dna_outcomes]
