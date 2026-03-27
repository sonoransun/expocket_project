"""2D editing map: sequence ruler + tool target-site bars + DICER contact track."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, Rectangle
from PySide6.QtCore import Signal

from viewer.config import NT_COLORS, EDITING_TOOL_COLORS, DOMAIN_COLORS

if TYPE_CHECKING:
    from viewer.data.schema import EnrichedVariantDataset
    from viewer.genome_editing.tools import TargetSite

# Rows in the map (from top):
# 0: sequence ruler
# 1: RNA track
# 2: DNA sense track
# 3: DICER pocket contact intensity track
# 4+: one row per tool

_TOOL_ORDER = ["cas9", "cas12a", "casclover", "talen", "zfn", "retron", "nicer"]
_TOOL_LABELS = {
    "cas9": "Cas9",
    "cas12a": "Cas12a",
    "casclover": "Cas-CLOVER",
    "talen": "TALEN",
    "zfn": "ZFN",
    "retron": "Retron/RLR",
    "nicer": "NICER",
}

_BG = "#1e1e1e"
_FG = "#cccccc"
_GRID = "#333333"


def plot_editing_map(
    fig: Figure,
    variant_id: str | None,
    rna_seq: str,
    sense_dna: str,
    sites_by_tool: dict[str, list["TargetSite"]],
    pocket_contact_vector: np.ndarray | None = None,
    highlighted_site: "TargetSite | None" = None,
) -> None:
    """Draw the full editing map for one variant."""
    fig.clear()
    fig.patch.set_facecolor(_BG)

    if not rna_seq:
        ax = fig.add_subplot(111)
        ax.set_facecolor(_BG)
        ax.text(0.5, 0.5, "Select a variant to see editing sites",
                ha="center", va="center", color=_FG, fontsize=9)
        return

    n_tools = len(_TOOL_ORDER)
    n_rows = 4 + n_tools  # ruler + RNA + DNA + contacts + N tools
    seq_len = len(rna_seq)

    # Build a GridSpec: rows have different heights
    row_heights = [0.4, 0.5, 0.5, 0.6] + [0.8] * n_tools
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(
        n_rows, 1,
        figure=fig,
        height_ratios=row_heights,
        hspace=0.1,
        left=0.12, right=0.98, top=0.93, bottom=0.06,
    )

    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_rows)]
    for ax in axes:
        ax.set_facecolor(_BG)
        ax.set_xlim(-0.5, seq_len - 0.5)
        for spine in ax.spines.values():
            spine.set_color(_GRID)
        ax.tick_params(colors=_FG, labelsize=5)

    # ── Row 0: sequence ruler ──────────────────────────────────────────────
    ax_ruler = axes[0]
    ax_ruler.set_ylim(0, 1)
    ax_ruler.set_yticks([])
    ax_ruler.set_xticks(range(0, seq_len, 5))
    ax_ruler.set_xticklabels([str(i) for i in range(0, seq_len, 5)], fontsize=5)
    ax_ruler.xaxis.tick_top()
    ax_ruler.xaxis.set_label_position("top")
    for x in range(0, seq_len, 5):
        ax_ruler.axvline(x, color=_GRID, lw=0.4, alpha=0.5)
    # Cleavage zone shading
    for xc, col in [(19, "#ca5cdd"), (20, "#5AC9A1"), (21, "#f94449"), (22, "#FFAE42")]:
        if xc < seq_len:
            ax_ruler.axvline(xc, color=col, lw=1.2, alpha=0.7, linestyle="--")
    title = f"Editing Site Map — {variant_id}" if variant_id else "Editing Site Map"
    fig.suptitle(title, color=_FG, fontsize=8, y=0.97)

    # ── Row 1: RNA sequence ────────────────────────────────────────────────
    _draw_sequence_track(axes[1], rna_seq, "RNA", use_rna_colors=True)

    # ── Row 2: DNA sense strand ────────────────────────────────────────────
    _draw_sequence_track(axes[2], sense_dna, "DNA", use_rna_colors=False)

    # ── Row 3: DICER pocket contact intensity ─────────────────────────────
    ax_contact = axes[3]
    ax_contact.set_yticks([])
    ax_contact.set_ylabel("DICER\ncontact", color=_FG, fontsize=5, rotation=0,
                          labelpad=40, va="center")
    if pocket_contact_vector is not None and len(pocket_contact_vector) > 0:
        cv = pocket_contact_vector[:seq_len]
        x = np.arange(len(cv))
        ax_contact.bar(x, cv, color="#5AC9A1", alpha=0.6, width=0.9)
        ax_contact.set_ylim(0, max(cv.max(), 0.01) * 1.2)
    else:
        ax_contact.text(seq_len / 2, 0.5, "No pocket data", ha="center",
                        va="center", color=_FG, fontsize=6, alpha=0.5)
        ax_contact.set_ylim(0, 1)

    # ── Rows 4+: tool site bars ────────────────────────────────────────────
    for i, tool_name in enumerate(_TOOL_ORDER):
        ax_tool = axes[4 + i]
        ax_tool.set_ylim(-0.1, 1.1)
        ax_tool.set_yticks([])
        label = _TOOL_LABELS.get(tool_name, tool_name)
        ax_tool.set_ylabel(label, color=_FG, fontsize=5, rotation=0,
                           labelpad=48, va="center")

        sites = sites_by_tool.get(tool_name, [])
        if not sites:
            ax_tool.text(seq_len / 2, 0.5, "no sites", ha="center",
                         va="center", color="#666", fontsize=5, style="italic")
            continue

        rgba = EDITING_TOOL_COLORS.get(tool_name, (0.7, 0.7, 0.7, 0.7))
        # Draw top-N sites as horizontal bars
        for site in sites[:10]:
            span_start = site.spacer_start
            span_end = site.spacer_end
            width = max(span_end - span_start, 1)
            color = list(rgba[:3]) + [0.6]
            rect = Rectangle(
                (span_start - 0.4, 0.2), width, 0.6,
                facecolor=color, edgecolor="none",
            )
            ax_tool.add_patch(rect)

            # PAM marker (small triangle)
            if site.pam_seq:
                pam_x = span_end if site.strand == "sense" else span_start
                ax_tool.annotate(
                    "▲", xy=(pam_x, 0.85), fontsize=4,
                    ha="center", color=list(rgba[:3]) + [1.0],
                )

            # Cut-site line
            if site.cut_positions[0] >= 0:
                cx = site.cut_positions[0]
                ax_tool.axvline(cx, color="white", lw=0.7, alpha=0.5, ymin=0.15, ymax=0.85)

            # Highlight selected site
            if (highlighted_site is not None and
                    highlighted_site.tool == tool_name and
                    highlighted_site.spacer_start == site.spacer_start and
                    highlighted_site.strand == site.strand):
                rect2 = Rectangle(
                    (span_start - 0.4, 0.1), width + 0.8, 0.8,
                    facecolor="none", edgecolor="white", lw=1.5,
                )
                ax_tool.add_patch(rect2)

    # Remove x-ticks from all rows except ruler
    for ax in axes[1:]:
        ax.set_xticks([])


def _draw_sequence_track(ax, seq: str, label: str, use_rna_colors: bool) -> None:
    """Draw a single-row sequence track with colored letter boxes."""
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel(label, color=_FG, fontsize=5, rotation=0, labelpad=28, va="center")
    for i, nt in enumerate(seq):
        key = nt.upper()
        if use_rna_colors:
            color = NT_COLORS.get(key, "#888888")
        else:
            # DNA: same but T instead of U
            color = NT_COLORS.get("T" if key == "T" else key, "#888888")
        ax.text(i, 0.5, nt.upper(), ha="center", va="center",
                fontsize=4, color=color, fontweight="bold")


class EditingMapWidget(FigureCanvasQTAgg):
    """Matplotlib canvas hosting the gene editing site map."""

    site_clicked = Signal(object)  # emits TargetSite

    def __init__(self, parent=None) -> None:
        self._fig = Figure(figsize=(8, 6), dpi=100)
        self._fig.set_facecolor(_BG)
        super().__init__(self._fig)
        self.setParent(parent)

        self._variant_id: str | None = None
        self._rna_seq: str = ""
        self._sense_dna: str = ""
        self._sites_by_tool: dict = {}
        self._pocket_contact: np.ndarray | None = None
        self._highlighted: "TargetSite | None" = None

        self.mpl_connect("button_press_event", self._on_click)

    def update_variant(
        self,
        variant_id: str,
        dataset: "EnrichedVariantDataset",
        sites_by_tool: dict[str, list["TargetSite"]],
    ) -> None:
        from viewer.genome_editing.sequence_utils import get_dual_strand
        variant = dataset.get_variant(variant_id)
        if variant is None:
            return
        self._variant_id = variant_id
        self._rna_seq = variant.pre_mirna_sequence
        sense, _ = get_dual_strand(self._rna_seq)
        self._sense_dna = sense
        self._sites_by_tool = sites_by_tool

        # Build contact intensity vector from pocket model
        self._pocket_contact = None
        if dataset.dicer_pocket is not None:
            seq_len = len(self._rna_seq)
            contact = np.zeros(seq_len, dtype=np.float64)
            for residue in dataset.dicer_pocket.residues:
                for pos in residue.rna_contact_positions:
                    if 0 <= pos < seq_len:
                        contact[pos] += residue.interaction_strength
            self._pocket_contact = contact

        self._redraw()

    def highlight_site(self, site: "TargetSite | None") -> None:
        self._highlighted = site
        self._redraw()

    def _redraw(self) -> None:
        plot_editing_map(
            self._fig,
            self._variant_id,
            self._rna_seq,
            self._sense_dna,
            self._sites_by_tool,
            self._pocket_contact,
            self._highlighted,
        )
        for ax in self._fig.get_axes():
            ax.set_facecolor(_BG)
            ax.tick_params(colors=_FG, labelsize=5)
            for spine in ax.spines.values():
                spine.set_color(_GRID)
        self.draw_idle()

    def _on_click(self, event) -> None:
        if event.inaxes is None:
            return
        # Identify which tool row was clicked
        axes = self._fig.get_axes()
        tool_axes = axes[4:]  # first 4 are ruler/RNA/DNA/contacts
        for i, ax in enumerate(tool_axes):
            if ax == event.inaxes:
                tool_name = _TOOL_ORDER[i] if i < len(_TOOL_ORDER) else None
                if tool_name is None:
                    return
                click_x = event.xdata
                sites = self._sites_by_tool.get(tool_name, [])
                best = None
                best_dist = float("inf")
                for site in sites[:10]:
                    mid = (site.spacer_start + site.spacer_end) / 2
                    d = abs(click_x - mid)
                    if d < best_dist:
                        best_dist = d
                        best = site
                if best is not None:
                    self.highlight_site(best)
                    self.site_clicked.emit(best)
                return
