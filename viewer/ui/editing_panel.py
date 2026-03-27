"""Gene editing dock panel — full pipeline: target sites, guide design, HDR, impact."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from viewer.genome_editing.tools import TOOLS, TOOL_NAMES, TargetSite, EditDesign
from viewer.genome_editing.sequence_utils import get_dual_strand

if TYPE_CHECKING:
    from viewer.data.schema import EnrichedVariantDataset
    from viewer.genome_editing.edit_engine import EditEngine
    from viewer.genome_editing.guide_designer import GuideDesigner
    from viewer.genome_editing.impact_predictor import EditImpactPredictor
    from viewer.genome_editing.off_target import OffTargetScorer
    from viewer.genome_editing.target_finder import TargetFinder

_DNA_BASES = ["A", "T", "G", "C"]


class EditingPanel(QWidget):
    """Dock panel for gene editing workflow: target → design → impact."""

    edit_designed = Signal(str, object)   # variant_id, EditDesign
    edit_applied = Signal(str, str)       # variant_id, modified_rna

    def __init__(
        self,
        dataset: "EnrichedVariantDataset",
        target_finder: "TargetFinder",
        guide_designer: "GuideDesigner",
        edit_engine: "EditEngine",
        impact_predictor: "EditImpactPredictor",
        off_target_scorer: "OffTargetScorer",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._dataset = dataset
        self._finder = target_finder
        self._guide_designer = guide_designer
        self._engine = edit_engine
        self._impact = impact_predictor
        self._ot_scorer = off_target_scorer

        self._current_variant: str | None = None
        self._sites: list[TargetSite] = []
        self._selected_site: TargetSite | None = None
        self._last_design: EditDesign | None = None
        self._pocket_only: bool = False

        layout = QVBoxLayout(self)
        layout.setSpacing(6)

        # ── Tool + sequence context ──────────────────────────────────────
        top_box = QGroupBox("Tool & Sequence")
        top_layout = QHBoxLayout(top_box)

        top_layout.addWidget(QLabel("Tool:"))
        self._tool_combo = QComboBox()
        for name in TOOL_NAMES:
            self._tool_combo.addItem(TOOLS[name].display_name, userData=name)
        self._tool_combo.currentIndexChanged.connect(self._on_tool_changed)
        top_layout.addWidget(self._tool_combo)

        top_layout.addWidget(QLabel("Seq:"))
        self._seq_group = QButtonGroup(self)
        self._dna_radio = QRadioButton("DNA")
        self._rna_radio = QRadioButton("RNA")
        self._dna_radio.setChecked(True)
        self._seq_group.addButton(self._dna_radio)
        self._seq_group.addButton(self._rna_radio)
        top_layout.addWidget(self._dna_radio)
        top_layout.addWidget(self._rna_radio)
        layout.addWidget(top_box)

        # ── Target sites ─────────────────────────────────────────────────
        sites_box = QGroupBox("Target Sites")
        sites_layout = QVBoxLayout(sites_box)

        self._sites_list = QListWidget()
        self._sites_list.setMaximumHeight(120)
        self._sites_list.currentRowChanged.connect(self._on_site_selected)
        sites_layout.addWidget(self._sites_list)

        btn_row = QHBoxLayout()
        rescan_btn = QPushButton("Rescan")
        rescan_btn.clicked.connect(self._rescan)
        btn_row.addWidget(rescan_btn)
        self._pocket_only_btn = QPushButton("Pocket Only")
        self._pocket_only_btn.setCheckable(True)
        self._pocket_only_btn.setChecked(False)
        self._pocket_only_btn.toggled.connect(self._on_pocket_only_toggled)
        btn_row.addWidget(self._pocket_only_btn)
        sites_layout.addLayout(btn_row)
        layout.addWidget(sites_box)

        # ── Sequence display ─────────────────────────────────────────────
        seq_box = QGroupBox("Sequence (selected context)")
        seq_layout = QVBoxLayout(seq_box)
        self._seq_display = QTextEdit()
        self._seq_display.setReadOnly(True)
        self._seq_display.setMaximumHeight(55)
        self._seq_display.setStyleSheet(
            "QTextEdit { font-family: monospace; font-size: 9px; "
            "background: #111; color: #aaa; border: 1px solid #333; }"
        )
        seq_layout.addWidget(self._seq_display)
        layout.addWidget(seq_box)

        # ── Edit design ──────────────────────────────────────────────────
        edit_box = QGroupBox("Edit Design")
        edit_layout = QVBoxLayout(edit_box)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Type:"))
        self._edit_type_combo = QComboBox()
        self._edit_type_combo.addItems(["snp", "insertion", "deletion", "hdr"])
        self._edit_type_combo.currentTextChanged.connect(self._on_edit_type_changed)
        row1.addWidget(self._edit_type_combo)
        row1.addWidget(QLabel("Pos:"))
        self._pos_spinner = QSpinBox()
        self._pos_spinner.setRange(0, 200)
        row1.addWidget(self._pos_spinner)
        edit_layout.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Seq/base:"))
        self._edit_seq = QLineEdit()
        self._edit_seq.setPlaceholderText("e.g. A or ACGU")
        row2.addWidget(self._edit_seq)
        design_btn = QPushButton("Design")
        design_btn.clicked.connect(self._run_design)
        row2.addWidget(design_btn)
        edit_layout.addLayout(row2)
        layout.addWidget(edit_box)

        # ── Results ──────────────────────────────────────────────────────
        results_box = QGroupBox("Results")
        results_layout = QVBoxLayout(results_box)

        def _lbl_row(label_text: str) -> tuple[QHBoxLayout, QLabel]:
            row = QHBoxLayout()
            row.addWidget(QLabel(label_text))
            val_lbl = QLabel("—")
            val_lbl.setStyleSheet("color: #5AC9A1; font-family: monospace; font-size: 8px;")
            row.addWidget(val_lbl)
            row.addStretch()
            results_layout.addLayout(row)
            return row, val_lbl

        _, self._guide_lbl = _lbl_row("Guide/Arms:")
        _, self._scaffold_lbl = _lbl_row("Scaffold/DR:")
        _, self._hdr_lbl = _lbl_row("HDR template:")
        _, self._retron_lbl = _lbl_row("Retron template:")
        _, self._ot_lbl = _lbl_row("Off-targets:")
        _, self._spec_lbl = _lbl_row("Specificity:")
        _, self._dc21_lbl = _lbl_row("DICER ΔDC21:")
        _, self._dc22_lbl = _lbl_row("DICER ΔDC22:")
        _, self._ratio_lbl = _lbl_row("Δratio (21-22):")
        _, self._struct_lbl = _lbl_row("Note:")
        _, self._cost_lbl = _lbl_row("Synthesis cost:")

        row_btns = QHBoxLayout()
        apply_btn = QPushButton("Apply to RNA Viewer")
        apply_btn.clicked.connect(self._apply_to_viewer)
        row_btns.addWidget(apply_btn)
        export_btn = QPushButton("Export JSON")
        export_btn.clicked.connect(self._export_json)
        row_btns.addWidget(export_btn)
        results_layout.addLayout(row_btns)
        layout.addWidget(results_box)
        layout.addStretch()

    # ── Public slots ────────────────────────────────────────────────────
    def set_variant(self, variant_id: str) -> None:
        self._current_variant = variant_id
        self._rescan()

    def select_site(self, site: "TargetSite") -> None:
        """Switch to site.tool, rescan, then select the matching row in the list."""
        combo_idx = self._tool_combo.findData(site.tool)
        if combo_idx >= 0 and combo_idx != self._tool_combo.currentIndex():
            self._tool_combo.blockSignals(True)
            self._tool_combo.setCurrentIndex(combo_idx)
            self._tool_combo.blockSignals(False)
        if self._current_variant is not None:
            self._rescan()
        for row in range(self._sites_list.count()):
            item = self._sites_list.item(row)
            if item is None:
                continue
            stored = item.data(256)
            if (stored is not None
                    and stored.spacer_start == site.spacer_start
                    and stored.strand == site.strand):
                self._sites_list.setCurrentRow(row)
                return

    # ── Internal ────────────────────────────────────────────────────────
    def _tool_name(self) -> str:
        idx = self._tool_combo.currentIndex()
        return self._tool_combo.itemData(idx) or TOOL_NAMES[0]

    def _on_tool_changed(self) -> None:
        if self._current_variant:
            self._rescan()

    def _on_pocket_only_toggled(self, checked: bool) -> None:
        self._pocket_only = checked
        self._populate_sites_list()

    def _on_edit_type_changed(self, edit_type: str) -> None:
        if edit_type == "snp":
            self._edit_seq.setPlaceholderText("single base, e.g. A")
        elif edit_type == "insertion":
            self._edit_seq.setPlaceholderText("bases to insert, e.g. ACGT")
        elif edit_type == "deletion":
            self._edit_seq.setPlaceholderText("bases to delete (number or sequence)")
        else:
            self._edit_seq.setPlaceholderText("new sequence for HDR replacement")

    def _rescan(self) -> None:
        if self._current_variant is None:
            return
        tool = self._tool_name()
        self._sites = self._finder.find_sites(self._current_variant, tool, self._dataset)
        self._populate_sites_list()

    def _populate_sites_list(self) -> None:
        self._sites_list.clear()
        pocket = getattr(self._dataset, "dicer_pocket", None)
        contact_positions: set[int] = set()
        if pocket is not None:
            for residue in pocket.residues:
                contact_positions.update(residue.rna_contact_positions)

        pocket_only = self._pocket_only
        visible_sites = []
        for site in self._sites[:20]:
            overlap = sum(
                1 for pos in range(site.spacer_start, site.spacer_end)
                if pos in contact_positions
            )
            if pocket_only and overlap == 0:
                continue
            visible_sites.append((site, overlap))

        for site, overlap in visible_sites:
            star = f"  \u2605{overlap}" if overlap > 0 else ""
            item = QListWidgetItem(f"{site.label()}{star}")
            item.setData(256, site)
            self._sites_list.addItem(item)

        if not visible_sites:
            msg = "(no pocket-overlap sites)" if pocket_only else "(no sites found)"
            self._sites_list.addItem(QListWidgetItem(msg))

    def _on_site_selected(self, row: int) -> None:
        if row < 0 or row >= len(self._sites):
            return
        self._selected_site = self._sites[row]
        # Update position spinner to cut site
        cut = self._selected_site.cut_positions[0]
        if cut >= 0:
            self._pos_spinner.setValue(cut)
        # Show sequence context
        self._update_seq_display()

    def _update_seq_display(self) -> None:
        if self._selected_site is None or self._current_variant is None:
            return
        variant = self._dataset.get_variant(self._current_variant)
        if variant is None:
            return
        rna = variant.pre_mirna_sequence
        sense, _ = get_dual_strand(rna)
        s = self._selected_site
        context_start = max(0, s.spacer_start - 3)
        context_end = min(len(sense), s.spacer_end + 3)
        if self._dna_radio.isChecked():
            display_seq = sense[context_start:context_end]
        else:
            display_seq = rna[context_start:context_end]
        label = f"pos {context_start}–{context_end}  |  spacer: {s.spacer_seq or 'N/A'}"
        self._seq_display.setPlainText(f"{label}\n{display_seq}")

    def _run_design(self) -> None:
        from viewer.genome_editing import guide_designer as gd
        if self._selected_site is None or self._current_variant is None:
            return
        site = self._selected_site
        variant = self._dataset.get_variant(self._current_variant)
        if variant is None:
            return

        original_rna = variant.pre_mirna_sequence
        sense_dna, _ = get_dual_strand(original_rna)
        edit_type = self._edit_type_combo.currentText()
        edit_seq = self._edit_seq.text().strip().upper() or "A"
        edit_pos = self._pos_spinner.value()
        tool_spec = TOOLS[site.tool]

        # --- Guide / binding domain sequences ---
        guide_seq = ""
        scaffold_str = ""
        guide2 = ""
        scaffold2_str = ""
        binding_arms = {}

        t = site.tool
        if t == "cas9":
            guide_seq, scaffold_str = gd.design_cas9_guide(site)
        elif t == "cas12a":
            scaffold_str, guide_seq = gd.design_cas12a_guide(site)
        elif t in ("casclover", "nicer"):
            guide_seq, scaffold_str, guide2, scaffold2_str = gd.design_casclover_pair(site)
        elif t == "talen":
            binding_arms = gd.design_talen_arms(site)
        elif t == "zfn":
            binding_arms = gd.design_zfn_fingers(site)
        elif t == "retron":
            guide_seq = ""

        # --- Edit simulation ---
        cut_pos = site.cut_positions[0] if site.cut_positions[0] >= 0 else edit_pos
        pam_pat = tool_spec.pam_pattern
        modified_rna, hdr_template = self._engine.apply_edit(
            original_rna, edit_type, edit_pos, edit_seq, cut_pos, pam_pat
        )
        retron_template = ""
        if t == "retron":
            retron_template = self._engine.retron_template(
                original_rna, edit_pos, edit_seq[:1] if edit_seq else "A"
            )
        nhej_outcomes = []
        if tool_spec.cut_type not in ("no_cut", "nick"):
            nhej_outcomes = self._engine.nhej_outcomes(original_rna, cut_pos)

        # --- Off-target scoring ---
        off_target_score, off_targets = 1.0, []
        if guide_seq:
            off_target_score, off_targets = self._ot_scorer.score(
                guide_seq, pam_pat,
                tool_spec.pam_position,
                tool_spec.spacer_len,
                self._dataset,
                exclude_variant_id=self._current_variant,
            )

        # --- DICER impact ---
        impact = self._impact.predict(self._current_variant, modified_rna)
        dicer_impact = impact.delta
        struct_note = impact.structural_note

        # --- Synthesis cost (rough: proportional to guide length × any mods) ---
        synthesis_cost = round(1.0 + len(guide_seq) * 0.02, 2)

        # Build EditDesign
        design = EditDesign(
            target_site=site,
            edit_type=edit_type,
            edit_sequence=edit_seq,
            edit_position=edit_pos,
            guide_sequence=guide_seq,
            scaffold=scaffold_str,
            guide2_sequence=guide2,
            scaffold2=scaffold2_str,
            binding_arms=binding_arms,
            hdr_template=hdr_template,
            retron_template=retron_template,
            nhej_outcomes=nhej_outcomes,
            off_target_score=off_target_score,
            off_targets=off_targets,
            predicted_modified_rna=modified_rna,
            dicer_impact=dicer_impact,
            synthesis_cost_factor=synthesis_cost,
            notes=[struct_note] if struct_note else [],
        )
        self._last_design = design
        self._populate_results(design)
        self.edit_designed.emit(self._current_variant, design)

    def _populate_results(self, design: EditDesign) -> None:
        guide = design.guide_sequence
        self._guide_lbl.setText(f"{guide[:30]}…" if len(guide) > 30 else (guide or "N/A"))
        sc = design.scaffold
        self._scaffold_lbl.setText(f"{sc[:25]}…" if len(sc) > 25 else (sc or "N/A"))
        hdr = design.hdr_template
        self._hdr_lbl.setText(f"{hdr[:25]}…" if len(hdr) > 25 else (hdr or "N/A"))
        ret = design.retron_template
        self._retron_lbl.setText(f"{ret[:25]}…" if len(ret) > 25 else (ret or "N/A"))
        n_ot = len(design.off_targets)
        self._ot_lbl.setText(f"{n_ot} in dataset")
        self._spec_lbl.setText(f"{design.off_target_score:.3f}")
        self._dc21_lbl.setText(f"{design.dicer_impact.get(21, 0.0):+.4f}")
        self._dc22_lbl.setText(f"{design.dicer_impact.get(22, 0.0):+.4f}")
        ratio = design.dicer_impact.get(21, 0.0) - design.dicer_impact.get(22, 0.0)
        self._ratio_lbl.setText(f"{ratio:+.4f}")
        self._struct_lbl.setText(design.notes[0] if design.notes else "—")
        self._cost_lbl.setText(f"{design.synthesis_cost_factor:.2f}×")

    def _apply_to_viewer(self) -> None:
        if self._last_design and self._current_variant:
            self.edit_applied.emit(
                self._current_variant,
                self._last_design.predicted_modified_rna,
            )

    def _export_json(self) -> None:
        if self._last_design is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Edit Design", "edit_design.json", "JSON (*.json)"
        )
        if not path:
            return
        site = self._last_design.target_site
        d = {
            "variant_id": site.variant_id,
            "tool": site.tool,
            "strand": site.strand,
            "spacer_start": site.spacer_start,
            "spacer_end": site.spacer_end,
            "spacer_seq": site.spacer_seq,
            "pam_seq": site.pam_seq,
            "cut_positions": list(self._last_design.target_site.cut_positions),
            "binding_score": site.binding_score,
            "edit_type": self._last_design.edit_type,
            "edit_sequence": self._last_design.edit_sequence,
            "edit_position": self._last_design.edit_position,
            "guide_sequence": self._last_design.guide_sequence,
            "scaffold": self._last_design.scaffold,
            "guide2_sequence": self._last_design.guide2_sequence,
            "binding_arms": self._last_design.binding_arms,
            "hdr_template": self._last_design.hdr_template,
            "retron_template": self._last_design.retron_template,
            "nhej_outcomes": self._last_design.nhej_outcomes,
            "off_target_score": self._last_design.off_target_score,
            "n_off_targets": len(self._last_design.off_targets),
            "dicer_impact": self._last_design.dicer_impact,
            "synthesis_cost_factor": self._last_design.synthesis_cost_factor,
            "notes": self._last_design.notes,
            "predicted_modified_rna": self._last_design.predicted_modified_rna,
        }
        with open(path, "w") as f:
            json.dump(d, f, indent=2)
