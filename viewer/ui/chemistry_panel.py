"""Chemistry lab dock panel for modification exploration."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from viewer.chemistry.modification_engine import ModificationEngine
    from viewer.chemistry.virtual_screen import ScreenResult, VirtualScreener
    from viewer.data.schema import EnrichedVariantDataset, ModificationState


class ChemistryPanel(QWidget):
    """Panel for applying chemical modifications and running virtual screens."""

    modification_applied = Signal(str, int, str)   # variant_id, position, mod_code
    modification_removed = Signal(str, int)         # variant_id, position
    screen_requested = Signal(str)                  # variant_id

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # Apply modification section
        apply_box = QGroupBox("Apply Modification")
        apply_layout = QVBoxLayout(apply_box)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Position:"))
        self._pos_spinner = QSpinBox()
        self._pos_spinner.setRange(0, 62)
        row1.addWidget(self._pos_spinner)
        row1.addWidget(QLabel("Mod:"))
        self._mod_combo = QComboBox()
        self._mod_combo.addItems(["2OMe", "LNA", "PSI", "m6A", "2F", "PS", "s4U", "m5C", "INO"])
        row1.addWidget(self._mod_combo)
        apply_layout.addLayout(row1)

        row2 = QHBoxLayout()
        self._apply_btn = QPushButton("Apply")
        self._apply_btn.clicked.connect(self._on_apply)
        row2.addWidget(self._apply_btn)
        self._remove_btn = QPushButton("Remove")
        self._remove_btn.clicked.connect(self._on_remove)
        row2.addWidget(self._remove_btn)
        self._clear_btn = QPushButton("Clear All")
        self._clear_btn.clicked.connect(self._on_clear)
        row2.addWidget(self._clear_btn)
        apply_layout.addLayout(row2)

        layout.addWidget(apply_box)

        # Current modifications table
        mod_box = QGroupBox("Current Modifications")
        mod_layout = QVBoxLayout(mod_box)
        self._mod_table = QTableWidget(0, 3)
        self._mod_table.setHorizontalHeaderLabels(["Position", "Nucleotide", "Modification"])
        self._mod_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._mod_table.setMaximumHeight(120)
        mod_layout.addWidget(self._mod_table)
        layout.addWidget(mod_box)

        # Virtual screen section
        screen_box = QGroupBox("Virtual Screen")
        screen_layout = QVBoxLayout(screen_box)
        self._screen_btn = QPushButton("Screen All Modifications")
        self._screen_btn.clicked.connect(self._on_screen)
        screen_layout.addWidget(self._screen_btn)

        self._screen_table = QTableWidget(0, 5)
        self._screen_table.setHorizontalHeaderLabels(
            ["Pos", "Mod", "ΔDC21", "ΔDC22", "ΔRatio"]
        )
        self._screen_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        screen_layout.addWidget(self._screen_table)
        layout.addWidget(screen_box)

        layout.addStretch()

        self._current_variant: str | None = None
        self._engine: ModificationEngine | None = None
        self._screener: VirtualScreener | None = None
        self._dataset: EnrichedVariantDataset | None = None

    def set_engines(
        self,
        engine: ModificationEngine,
        screener: VirtualScreener,
        dataset: EnrichedVariantDataset,
    ) -> None:
        self._engine = engine
        self._screener = screener
        self._dataset = dataset

    def set_variant(self, variant_id: str) -> None:
        self._current_variant = variant_id
        self._update_mod_table()

    def _on_apply(self) -> None:
        if self._current_variant and self._engine:
            pos = self._pos_spinner.value()
            mod = self._mod_combo.currentText()
            try:
                self._engine.apply_modification(self._current_variant, pos, mod)
                self._update_mod_table()
                self.modification_applied.emit(self._current_variant, pos, mod)
            except ValueError:
                pass

    def _on_remove(self) -> None:
        if self._current_variant and self._engine:
            pos = self._pos_spinner.value()
            self._engine.remove_modification(self._current_variant, pos)
            self._update_mod_table()
            self.modification_removed.emit(self._current_variant, pos)

    def _on_clear(self) -> None:
        if self._current_variant and self._engine:
            self._engine.clear_modifications(self._current_variant)
            self._update_mod_table()

    def _on_screen(self) -> None:
        if self._current_variant and self._screener:
            results = self._screener.rank_by_dc_ratio_shift(
                self._current_variant, top_n=15
            )
            self._populate_screen_table(results)

    def _update_mod_table(self) -> None:
        self._mod_table.setRowCount(0)
        if not self._current_variant or not self._dataset:
            return
        state = self._dataset.get_modification_state(self._current_variant)
        variant = self._dataset.get_variant(self._current_variant)
        if not variant:
            return
        for pos, mod_code in sorted(state.modifications.items()):
            row = self._mod_table.rowCount()
            self._mod_table.insertRow(row)
            nt = variant.pre_mirna_sequence[pos] if pos < len(variant.pre_mirna_sequence) else "?"
            self._mod_table.setItem(row, 0, QTableWidgetItem(str(pos)))
            self._mod_table.setItem(row, 1, QTableWidgetItem(nt))
            self._mod_table.setItem(row, 2, QTableWidgetItem(mod_code))

    def _populate_screen_table(self, results: list[ScreenResult]) -> None:
        self._screen_table.setRowCount(0)
        for r in results:
            row = self._screen_table.rowCount()
            self._screen_table.insertRow(row)
            self._screen_table.setItem(row, 0, QTableWidgetItem(str(r.position)))
            self._screen_table.setItem(row, 1, QTableWidgetItem(r.modification))
            self._screen_table.setItem(row, 2, QTableWidgetItem(f"{r.delta_dc21:+.4f}"))
            self._screen_table.setItem(row, 3, QTableWidgetItem(f"{r.delta_dc22:+.4f}"))
            self._screen_table.setItem(row, 4, QTableWidgetItem(f"{r.delta_ratio:+.4f}"))
