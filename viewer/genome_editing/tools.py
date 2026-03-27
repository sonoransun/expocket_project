"""Canonical gene editing tool specs and shared data structures."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ToolSpec:
    """Immutable specification for one gene editing tool."""
    name: str
    display_name: str
    pam_pattern: str | None       # IUPAC PAM; None for TALEN/ZFN/Retron
    pam_position: str             # "3prime" | "5prime" | "none"
    spacer_len: int               # guide/spacer length in nt (0 for TALEN/ZFN)
    cut_type: str                 # "blunt" | "staggered" | "nick" | "paired_nick" | "no_cut"
    cut_offset: int               # nt from 5' end of spacer to cut (sense strand)
    overhang: int                 # 0=blunt, >0 = 5' overhang nt, <0 = 3' overhang
    color: tuple[float, float, float, float]  # RGBA


TOOLS: dict[str, ToolSpec] = {
    "cas9": ToolSpec(
        name="cas9",
        display_name="CRISPR-Cas9 (SpCas9)",
        pam_pattern="NGG",
        pam_position="3prime",
        spacer_len=20,
        cut_type="blunt",
        cut_offset=17,   # cut between pos 17-18 from spacer 5' end (3 bp from PAM)
        overhang=0,
        color=(0.2, 0.6, 1.0, 0.85),
    ),
    "cas12a": ToolSpec(
        name="cas12a",
        display_name="Cas12a (Cpf1)",
        pam_pattern="TTTN",
        pam_position="5prime",
        spacer_len=21,
        cut_type="staggered",
        cut_offset=18,   # non-template strand cut 18 nt from PAM; template 23 → 5-nt 5' overhang
        overhang=5,
        color=(0.4, 0.9, 0.4, 0.85),
    ),
    "casclover": ToolSpec(
        name="casclover",
        display_name="Cas-CLOVER (paired Cas9n)",
        pam_pattern="NGG",
        pam_position="3prime",
        spacer_len=20,
        cut_type="paired_nick",
        cut_offset=17,
        overhang=25,     # window between the two guide pairs (25-50 bp)
        color=(1.0, 0.7, 0.1, 0.85),
    ),
    "talen": ToolSpec(
        name="talen",
        display_name="TALEN",
        pam_pattern=None,
        pam_position="none",
        spacer_len=0,
        cut_type="staggered",
        cut_offset=0,
        overhang=-1,     # 1-nt 3' overhang from FokI in spacer
        color=(0.9, 0.3, 0.3, 0.85),
    ),
    "zfn": ToolSpec(
        name="zfn",
        display_name="ZFN (Zinc Finger Nuclease)",
        pam_pattern=None,
        pam_position="none",
        spacer_len=0,
        cut_type="staggered",
        cut_offset=0,
        overhang=-4,     # 4-nt 3' overhang from FokI
        color=(0.8, 0.4, 0.9, 0.85),
    ),
    "retron": ToolSpec(
        name="retron",
        display_name="Retron/RLR (EC86-style)",
        pam_pattern=None,
        pam_position="none",
        spacer_len=0,
        cut_type="no_cut",
        cut_offset=0,
        overhang=0,
        color=(0.3, 0.9, 0.8, 0.85),
    ),
    "nicer": ToolSpec(
        name="nicer",
        display_name="NICER (paired D10A nicks)",
        pam_pattern="NGG",
        pam_position="3prime",
        spacer_len=20,
        cut_type="paired_nick",
        cut_offset=17,
        overhang=10,     # window between paired nicks: 5-25 bp
        color=(1.0, 0.5, 0.2, 0.85),
    ),
}

TOOL_NAMES: list[str] = list(TOOLS.keys())
TOOL_DISPLAY_NAMES: list[str] = [t.display_name for t in TOOLS.values()]


@dataclass
class TargetSite:
    """One valid target site for a gene editing tool on a variant sequence."""
    tool: str
    variant_id: str
    strand: str               # "sense" | "antisense"
    spacer_start: int         # sense-strand coordinate, 0-based
    spacer_end: int
    spacer_seq: str           # DNA spacer (20-21 nt, or "" for TALEN/ZFN)
    pam_seq: str              # actual PAM bases, or "" for TALEN/ZFN/Retron
    cut_positions: tuple[int, int]  # (cut_sense_strand, cut_antisense_strand) sense coords
    gc_content: float         # spacer GC fraction
    seed_score: float         # 0-1 guide seed-region quality
    binding_score: float      # composite quality score 0-1
    # TALEN/ZFN extras
    binding_arms: dict = field(default_factory=dict)
    # Paired sites (Cas-CLOVER / NICER)
    paired_site: "TargetSite | None" = field(default=None, repr=False)

    def label(self) -> str:
        strand_sym = "+" if self.strand == "sense" else "-"
        tool_short = {"cas9": "Cas9", "cas12a": "Cas12a", "casclover": "CasC",
                      "talen": "TALEN", "zfn": "ZFN", "retron": "Ret", "nicer": "NICER"}
        return f"{tool_short.get(self.tool, self.tool)}  {strand_sym}{self.spacer_start}–{self.spacer_end}  score={self.binding_score:.2f}"


@dataclass
class EditDesign:
    """Complete gene editing design for a single target site."""
    target_site: TargetSite
    edit_type: str            # "snp" | "insertion" | "deletion" | "hdr"
    edit_sequence: str        # new base(s), or "" for deletion
    edit_position: int        # sense-strand DNA coordinate of the edit
    # Guide / binding domain sequences
    guide_sequence: str       # 20-21 nt spacer (DNA), or "" for TALEN/ZFN
    scaffold: str             # sgRNA scaffold or crRNA direct repeat, or ""
    guide2_sequence: str      # second guide for paired tools, or ""
    scaffold2: str            # second scaffold
    # TALEN/ZFN arm info stored in target_site.binding_arms
    # Templates
    hdr_template: str         # ~100 nt ssODN for HDR
    retron_template: str      # EC86 msd template, or ""
    nhej_outcomes: list[str]  # top predicted NHEJ indels (for DSB tools)
    # Scores
    off_target_score: float   # 0-1 (1 = perfectly specific within dataset)
    off_targets: list[dict]   # [{"variant_id":..,"mismatches":..,"site":..}]
    # Predicted biological outcome
    predicted_modified_rna: str
    dicer_impact: dict        # {20: delta, 21: delta, 22: delta, 23: delta}
    synthesis_cost_factor: float
    notes: list[str]          # warnings, incompatibilities
