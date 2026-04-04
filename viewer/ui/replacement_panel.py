"""Single and double replacement chemistry workflow panel."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from PySide6.QtCore import Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
        from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
    from viewer.chemistry.double_screen import DoubleReplacementScreener, DoubleScreenResult
    from viewer.chemistry.modification_engine import ModificationEngine
    from viewer.data.schema import EnrichedVariantDataset

from viewer.config import SYNERGY_COLORS

_log = logging.getLogger(__name__)

_MOD_OPTIONS = ["A", "U", "G", "C", "2OMe", "LNA", "PSI", "m6A", "2F", "PS", "s4U", "m5C", "INO"]
_BASE_CHANGES = {"A", "U", "G", "C"}


class ReplacementPanel(QWidget):
    """Panel for single and double nucleotide/modification replacement."""

    replacement_applied = Signal(str, dict)   # variant_id, {pos: mod_code}
    comparison_ready = Signal(dict)           # comparison data for 2D plot
    add_to_batch_requested = Signal(str, str) # variant_id, mode ("double")

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setSpacing(6)

        self._tabs = QTabWidget()

        # --- Single replacement tab ---
        single_widget = QWidget()
        single_layout = QVBoxLayout(single_widget)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Position:"))
        self._single_pos = QSpinBox()
        self._single_pos.setRange(0, 200)
        row1.addWidget(self._single_pos)
        row1.addWidget(QLabel("Current:"))
        self._single_current_nt = QLabel("—")
        row1.addWidget(self._single_current_nt)
        row1.addWidget(QLabel("Replace with:"))
        self._single_replacement = QComboBox()
        self._single_replacement.addItems(_MOD_OPTIONS)
        row1.addWidget(self._single_replacement)
        single_layout.addLayout(row1)

        row2 = QHBoxLayout()
        self._single_preview_btn = QPushButton("Preview")
        self._single_preview_btn.clicked.connect(self._on_single_preview)
        row2.addWidget(self._single_preview_btn)
        self._single_apply_btn = QPushButton("Apply")
        self._single_apply_btn.clicked.connect(self._on_single_apply)
        row2.addWidget(self._single_apply_btn)
        single_layout.addLayout(row2)

        self._single_result_table = QTableWidget(4, 3)
        self._single_result_table.setHorizontalHeaderLabels(["Site", "Original", "Modified"])
        self._single_result_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._single_result_table.setMaximumHeight(130)
        single_layout.addWidget(self._single_result_table)
        single_layout.addStretch()

        self._tabs.addTab(single_widget, "Single")

        # --- Double replacement tab ---
        double_widget = QWidget()
        double_layout = QVBoxLayout(double_widget)

        # Position 1
        r1 = QHBoxLayout()
        r1.addWidget(QLabel("Pos 1:"))
        self._double_pos1 = QSpinBox()
        self._double_pos1.setRange(0, 200)
        r1.addWidget(self._double_pos1)
        r1.addWidget(QLabel("Mod 1:"))
        self._double_mod1 = QComboBox()
        self._double_mod1.addItems(_MOD_OPTIONS)
        r1.addWidget(self._double_mod1)
        double_layout.addLayout(r1)

        # Position 2
        r2 = QHBoxLayout()
        r2.addWidget(QLabel("Pos 2:"))
        self._double_pos2 = QSpinBox()
        self._double_pos2.setRange(0, 200)
        self._double_pos2.setValue(5)
        r2.addWidget(self._double_pos2)
        r2.addWidget(QLabel("Mod 2:"))
        self._double_mod2 = QComboBox()
        self._double_mod2.addItems(_MOD_OPTIONS)
        r2.addWidget(self._double_mod2)
        double_layout.addLayout(r2)

        # Screen controls
        screen_row = QHBoxLayout()
        self._screen_btn = QPushButton("Screen Combinations")
        self._screen_btn.clicked.connect(self._on_screen)
        screen_row.addWidget(self._screen_btn)
        self._clv_zone_check = QCheckBox("Cleavage zone only")
        self._clv_zone_check.setChecked(True)
        screen_row.addWidget(self._clv_zone_check)
        self._batch_btn = QPushButton("Add to Batch")
        self._batch_btn.clicked.connect(self._on_add_to_batch)
        screen_row.addWidget(self._batch_btn)
        double_layout.addLayout(screen_row)

        # Results table
        self._double_table = QTableWidget(0, 7)
        self._double_table.setHorizontalHeaderLabels(
            ["Pos1", "Mod1", "Pos2", "Mod2", "ΔRatio", "Synergy", "Yield%"]
        )
        self._double_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._double_table.cellClicked.connect(self._on_double_row_clicked)
        double_layout.addWidget(self._double_table)

        apply_row = QHBoxLayout()
        self._double_preview_btn = QPushButton("Preview Selected")
        self._double_preview_btn.clicked.connect(self._on_double_preview)
        apply_row.addWidget(self._double_preview_btn)
        self._double_apply_btn = QPushButton("Apply Selected")
        self._double_apply_btn.clicked.connect(self._on_double_apply)
        apply_row.addWidget(self._double_apply_btn)
        double_layout.addLayout(apply_row)

        self._tabs.addTab(double_widget, "Double")
        layout.addWidget(self._tabs)

        self._current_variant: str | None = None
        self._dataset: EnrichedVariantDataset | None = None
        self._engine: ModificationEngine | None = None
        self._predictor: CleavageSitePredictor | None = None
        self._double_screener: DoubleReplacementScreener | None = None
        self._double_results: list[DoubleScreenResult] = []
        self._selected_double_row: int = -1

    def set_engines(
        self,
        engine: ModificationEngine,
        predictor: CleavageSitePredictor,
        double_screener: DoubleReplacementScreener,
        dataset: EnrichedVariantDataset,
    ) -> None:
        self._engine = engine
        self._predictor = predictor
        self._double_screener = double_screener
        self._dataset = dataset

    def set_variant(self, variant_id: str) -> None:
        self._current_variant = variant_id
        # Update current nucleotide display and spinner ranges
        if self._dataset:
            vi = self._dataset.get_variant(variant_id)
            if vi:
                max_pos = max(0, len(vi.pre_mirna_sequence) - 1)
                self._single_pos.setRange(0, max_pos)
                self._double_pos1.setRange(0, max_pos)
                self._double_pos2.setRange(0, max_pos)
                pos = self._single_pos.value()
                if pos < len(vi.pre_mirna_sequence):
                    self._single_current_nt.setText(vi.pre_mirna_sequence[pos])

    def _get_single_mods(self) -> dict[int, str]:
        pos = self._single_pos.value()
        replacement = self._single_replacement.currentText()
        if replacement in _BASE_CHANGES:
            return {pos: f"BASE:{replacement}"}
        return {pos: replacement}

    def _on_single_preview(self) -> None:
        if not self._current_variant or not self._predictor:
            return
        mods = self._get_single_mods()
        if not mods:
            return
        pos, code = next(iter(mods.items()))
        if code.startswith("BASE:"):
            new_base = code.split(":")[1]
            shifts = self._predictor.predict_base_change(
                self._current_variant, pos, new_base
            )
        else:
            shifts = self._predictor.predict_shift(self._current_variant, mods)
        self._single_result_table.setRowCount(4)
        for i, site in enumerate([20, 21, 22, 23]):
            self._single_result_table.setItem(i, 0, QTableWidgetItem(f"DC{site}"))
            self._single_result_table.setItem(i, 1, QTableWidgetItem("0.000"))
            self._single_result_table.setItem(i, 2, QTableWidgetItem(f"{shifts.get(site, 0.0):+.4f}"))

        # Emit comparison for 2D plot
        self.comparison_ready.emit({
            "original": {s: 0.0 for s in [20, 21, 22, 23]},
            "single_1": shifts,
        })

    def _on_single_apply(self) -> None:
        if not self._current_variant or not self._engine:
            return
        mods = self._get_single_mods()
        if mods:
            for pos, mod in mods.items():
                if mod.startswith("BASE:"):
                    _log.info("Base changes are preview-only; use Gene Editing panel to apply sequence mutations")
                    return
                try:
                    self._engine.apply_modification(self._current_variant, pos, mod)
                except ValueError as exc:
                    _log.info("Modification rejected at pos %d: %s", pos, exc)
                    return
            self.replacement_applied.emit(self._current_variant, mods)

    def _on_add_to_batch(self) -> None:
        if self._current_variant:
            self.add_to_batch_requested.emit(self._current_variant, "double")

    def _on_screen(self) -> None:
        if not self._current_variant or not self._double_screener:
            return
        positions = None
        if not self._clv_zone_check.isChecked():
            vi = self._dataset.get_variant(self._current_variant) if self._dataset else None
            seq_len = len(vi.pre_mirna_sequence) if vi else 63
            positions = list(range(seq_len))

        self._double_results = self._double_screener.screen_double(
            self._current_variant, positions=positions, top_n=30
        )
        self._populate_double_table()

    def _populate_double_table(self) -> None:
        self._double_table.setRowCount(0)
        for r in self._double_results:
            row = self._double_table.rowCount()
            self._double_table.insertRow(row)
            self._double_table.setItem(row, 0, QTableWidgetItem(str(r.position_1)))
            self._double_table.setItem(row, 1, QTableWidgetItem(r.mod_1))
            self._double_table.setItem(row, 2, QTableWidgetItem(str(r.position_2)))
            self._double_table.setItem(row, 3, QTableWidgetItem(r.mod_2))
            self._double_table.setItem(row, 4, QTableWidgetItem(f"{r.delta_ratio:+.4f}"))

            # Color synergy
            syn_item = QTableWidgetItem(f"{r.synergy:+.4f}")
            if r.synergy > 0.001:
                syn_item.setForeground(QColor(SYNERGY_COLORS["cooperative"]))
            elif r.synergy < -0.001:
                syn_item.setForeground(QColor(SYNERGY_COLORS["antagonistic"]))
            else:
                syn_item.setForeground(QColor(SYNERGY_COLORS["neutral"]))
            self._double_table.setItem(row, 5, syn_item)

            self._double_table.setItem(row, 6, QTableWidgetItem(f"{r.synthesis_yield * 100:.1f}"))

    def _on_double_row_clicked(self, row: int, col: int) -> None:
        self._selected_double_row = row

    def _on_double_preview(self) -> None:
        if self._selected_double_row < 0 or self._selected_double_row >= len(self._double_results):
            return
        if not self._current_variant or not self._double_screener:
            return
        r = self._double_results[self._selected_double_row]
        cmp = self._double_screener.compare_single_vs_double(
            self._current_variant, r.position_1, r.mod_1, r.position_2, r.mod_2
        )
        self.comparison_ready.emit(cmp)

    def _on_double_apply(self) -> None:
        if self._selected_double_row < 0 or self._selected_double_row >= len(self._double_results):
            return
        if not self._current_variant or not self._engine:
            return
        r = self._double_results[self._selected_double_row]
        mods = {r.position_1: r.mod_1, r.position_2: r.mod_2}
        for pos, mod in mods.items():
            try:
                self._engine.apply_modification(self._current_variant, pos, mod)
            except ValueError as exc:
                _log.info("Modification rejected at pos %d: %s", pos, exc)
                return
        self.replacement_applied.emit(self._current_variant, mods)
