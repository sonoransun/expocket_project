"""Synthesis planning panel showing step-by-step oligonucleotide synthesis."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QComboBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from viewer.chemistry.synthesis_pathway import SynthesisPlanner
    from viewer.data.schema import EnrichedVariantDataset, SynthesisPlan


class SynthesisPanel(QWidget):
    """Panel showing step-by-step synthesis plan with yield and cost info."""

    synthesis_updated = Signal(str)  # variant_id

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setSpacing(6)

        # Header
        header = QHBoxLayout()
        self._variant_label = QLabel("No variant selected")
        self._variant_label.setStyleSheet("font-weight: bold; font-size: 12px;")
        header.addWidget(self._variant_label)
        header.addStretch()
        header.addWidget(QLabel("Scale:"))
        self._scale_combo = QComboBox()
        self._scale_combo.addItems(["50 nmol", "200 nmol", "1000 nmol"])
        self._scale_combo.setCurrentIndex(1)
        self._scale_combo.currentIndexChanged.connect(self._on_scale_changed)
        header.addWidget(self._scale_combo)
        layout.addLayout(header)

        # Summary bar
        summary_box = QGroupBox("Summary")
        summary_layout = QHBoxLayout(summary_box)
        self._yield_label = QLabel("Yield: —")
        self._cost_label = QLabel("Cost: —")
        self._length_label = QLabel("Length: —")
        self._mods_label = QLabel("Mods: —")
        for lbl in (self._yield_label, self._cost_label, self._length_label, self._mods_label):
            lbl.setStyleSheet("font-size: 11px;")
            summary_layout.addWidget(lbl)
        layout.addWidget(summary_box)

        # Step table
        self._step_table = QTableWidget(0, 9)
        self._step_table.setHorizontalHeaderLabels(
            ["Step", "Pos", "Base", "Mod", "Monomer", "Eff%", "Yield%", "Cost", "Notes"]
        )
        self._step_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._step_table.setMaximumHeight(250)
        layout.addWidget(self._step_table)

        # Warnings
        self._warnings_box = QTextEdit()
        self._warnings_box.setReadOnly(True)
        self._warnings_box.setMaximumHeight(60)
        self._warnings_box.setStyleSheet("color: #f44; font-size: 10px; background: #1e1e1e;")
        layout.addWidget(self._warnings_box)

        layout.addStretch()

        self._planner: SynthesisPlanner | None = None
        self._dataset: EnrichedVariantDataset | None = None
        self._current_variant: str | None = None
        self._current_plan: SynthesisPlan | None = None

    def set_planner(
        self, planner: SynthesisPlanner, dataset: EnrichedVariantDataset
    ) -> None:
        self._planner = planner
        self._dataset = dataset

    def set_variant(self, variant_id: str) -> None:
        self._current_variant = variant_id
        self._variant_label.setText(f"Synthesis: {variant_id}")
        self._refresh_plan()

    def refresh(self) -> None:
        """Recalculate synthesis plan for current variant."""
        self._refresh_plan()

    @property
    def current_plan(self) -> SynthesisPlan | None:
        return self._current_plan

    def _on_scale_changed(self) -> None:
        self._refresh_plan()

    def _get_scale(self) -> int:
        text = self._scale_combo.currentText()
        return int(text.split()[0])

    def _refresh_plan(self) -> None:
        if not self._planner or not self._current_variant or not self._dataset:
            return

        state = self._dataset.get_modification_state(self._current_variant)
        plan = self._planner.plan_synthesis(
            self._current_variant,
            dict(state.modifications),
            scale_nmol=self._get_scale(),
        )
        self._current_plan = plan
        self._update_display(plan)
        self.synthesis_updated.emit(self._current_variant)

    def _update_display(self, plan: SynthesisPlan) -> None:
        # Summary
        self._yield_label.setText(f"Yield: {plan.total_yield * 100:.1f}%")
        self._cost_label.setText(f"Cost: {plan.total_cost_factor:.0f}×")
        self._length_label.setText(f"Length: {len(plan.sequence)} nt")
        self._mods_label.setText(f"Mods: {len(plan.modifications)}")

        # Color yield label
        if plan.total_yield > 0.5:
            self._yield_label.setStyleSheet("color: #5AC9A1; font-size: 11px;")
        elif plan.total_yield > 0.2:
            self._yield_label.setStyleSheet("color: #FFAE42; font-size: 11px;")
        else:
            self._yield_label.setStyleSheet("color: #f44; font-size: 11px;")

        # Step table
        self._step_table.setRowCount(0)
        for i, step in enumerate(plan.steps):
            row = self._step_table.rowCount()
            self._step_table.insertRow(row)

            items = [
                str(i + 1),
                str(step.position),
                step.nucleotide,
                step.modification or "—",
                step.monomer_name,
                f"{step.coupling_efficiency * 100:.1f}",
                f"{step.cumulative_yield * 100:.1f}",
                f"{step.cost_factor:.1f}",
                step.notes,
            ]
            for col, text in enumerate(items):
                item = QTableWidgetItem(text)
                self._step_table.setItem(row, col, item)

            # Row coloring based on coupling efficiency
            if step.coupling_efficiency >= 0.99:
                bg = QColor(40, 80, 40)
            elif step.coupling_efficiency >= 0.98:
                bg = QColor(80, 80, 30)
            else:
                bg = QColor(80, 30, 30)

            if step.modification:
                bg = QColor(50, 50, 80)

            for col in range(9):
                item = self._step_table.item(row, col)
                if item:
                    item.setBackground(bg)

        # Warnings
        if plan.incompatibilities:
            self._warnings_box.setText("\n".join(plan.incompatibilities))
        else:
            self._warnings_box.setText("No compatibility issues detected.")
            self._warnings_box.setStyleSheet("color: #5AC9A1; font-size: 10px; background: #1e1e1e;")
