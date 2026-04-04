"""Batch screening panel: queue variants, run in background, view results live."""

from __future__ import annotations

import csv
from datetime import datetime
from typing import TYPE_CHECKING

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from viewer.chemistry.batch_screener import BatchJob, BatchResult
    from viewer.data.schema import EnrichedVariantDataset

_GROUPS = ["A", "T", "G", "C"]


class BatchPanel(QWidget):
    """Dock panel for batch multi-variant screening.

    Users queue variants individually (via "Add to Batch" from other panels)
    or in bulk ("Add Group X" buttons), then click "Run Batch".  Results
    arrive opportunistically as each variant completes.
    """

    run_requested = Signal(list)    # list[BatchJob] — picked up by MainWindow
    cancel_requested = Signal()
    jump_to_variant = Signal(str)   # variant_id → controller.select_variant

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setSpacing(6)

        # ── Job Queue ────────────────────────────────────────────────────────
        queue_box = QGroupBox("Job Queue")
        queue_layout = QVBoxLayout(queue_box)

        self._queue_table = QTableWidget(0, 3)
        self._queue_table.setHorizontalHeaderLabels(["Variant", "Mode", "Status"])
        self._queue_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self._queue_table.setMaximumHeight(130)
        queue_layout.addWidget(self._queue_table)

        # "Add Group X" buttons + mode selector
        add_row = QHBoxLayout()
        add_row.addWidget(QLabel("Add group:"))
        for grp in _GROUPS:
            btn = QPushButton(grp)
            btn.setMaximumWidth(36)
            btn.clicked.connect(lambda checked, g=grp: self._add_group(g))
            add_row.addWidget(btn)
        add_row.addWidget(QLabel("Mode:"))
        self._mode_combo = QComboBox()
        self._mode_combo.addItems(["single", "double"])
        add_row.addWidget(self._mode_combo)
        add_row.addWidget(QLabel("Top N:"))
        self._topn_spin = QSpinBox()
        self._topn_spin.setRange(5, 100)
        self._topn_spin.setValue(20)
        add_row.addWidget(self._topn_spin)
        add_row.addStretch()
        queue_layout.addLayout(add_row)

        # Run / Cancel / progress row
        ctrl_row = QHBoxLayout()
        self._run_btn = QPushButton("Run Batch")
        self._run_btn.clicked.connect(self._on_run)
        ctrl_row.addWidget(self._run_btn)
        self._cancel_btn = QPushButton("Cancel")
        self._cancel_btn.setEnabled(False)
        self._cancel_btn.clicked.connect(self.cancel_requested)
        ctrl_row.addWidget(self._cancel_btn)
        self._progress_label = QLabel("0 / 0")
        ctrl_row.addWidget(self._progress_label)
        ctrl_row.addStretch()
        queue_layout.addLayout(ctrl_row)

        layout.addWidget(queue_box)

        # ── Progress bar ─────────────────────────────────────────────────────
        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 100)
        self._progress_bar.setValue(0)
        layout.addWidget(self._progress_bar)

        # ── Live Results ──────────────────────────────────────────────────────
        results_box = QGroupBox("Results")
        results_layout = QVBoxLayout(results_box)

        self._results_table = QTableWidget(0, 6)
        self._results_table.setHorizontalHeaderLabels(
            ["Variant", "Pos", "Mod", "ΔDC21", "ΔDC22", "ΔRatio"]
        )
        self._results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        results_layout.addWidget(self._results_table)

        jump_row = QHBoxLayout()
        self._jump_btn = QPushButton("Jump to Variant")
        self._jump_btn.setEnabled(False)
        self._jump_btn.clicked.connect(self._on_jump)
        jump_row.addWidget(self._jump_btn)
        self._export_btn = QPushButton("Export CSV")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._export_csv)
        jump_row.addWidget(self._export_btn)
        jump_row.addStretch()
        results_layout.addLayout(jump_row)
        self._results_table.cellClicked.connect(
            lambda r, c: self._jump_btn.setEnabled(True)
        )

        layout.addWidget(results_box)

        # ── Aggregated Summary ───────────────────────────────────────────────
        summary_box = QGroupBox("Best Hit per Variant")
        summary_layout = QVBoxLayout(summary_box)
        self._summary_table = QTableWidget(0, 4)
        self._summary_table.setHorizontalHeaderLabels(
            ["Variant", "Best Mod", "ΔRatio", "Mode"]
        )
        self._summary_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self._summary_table.setMaximumHeight(150)
        self._summary_table.cellClicked.connect(
            lambda r, c: self._jump_btn.setEnabled(True)
        )
        summary_layout.addWidget(self._summary_table)
        layout.addWidget(summary_box)

        # Internal state
        self._pending_jobs: list[BatchJob] = []
        self._dataset: EnrichedVariantDataset | None = None
        # {variant_id: row_index} for queue table lookups
        self._queue_row: dict[str, int] = {}
        # {variant_id: best_delta_ratio} for summary upsert
        self._summary_best: dict[str, float] = {}
        # {variant_id: summary_row} for upsert
        self._summary_row: dict[str, int] = {}

    def set_dataset(self, dataset: EnrichedVariantDataset) -> None:
        self._dataset = dataset

    # ── Public API ────────────────────────────────────────────────────────────

    def add_job(self, variant_id: str, mode: str) -> None:
        """Add a single variant to the queue; silently skips duplicates."""
        from viewer.chemistry.batch_screener import BatchJob

        key = (variant_id, mode)
        if any((j.variant_id, j.mode) == key for j in self._pending_jobs):
            return

        job = BatchJob(
            variant_id=variant_id,
            mode=mode,
            top_n=self._topn_spin.value(),
        )
        self._pending_jobs.append(job)

        row = self._queue_table.rowCount()
        self._queue_table.insertRow(row)
        self._queue_table.setItem(row, 0, QTableWidgetItem(variant_id))
        self._queue_table.setItem(row, 1, QTableWidgetItem(mode))
        self._queue_table.setItem(row, 2, QTableWidgetItem("queued"))
        self._queue_row[variant_id] = row

    # ── Slots called by BatchWorker signals (main thread via QueuedConnection) ─

    def _on_result(self, result: BatchResult) -> None:
        # Update queue table status
        row_idx = self._queue_row.get(result.variant_id)
        if row_idx is not None:
            status = "error" if result.error else "done"
            self._queue_table.setItem(row_idx, 2, QTableWidgetItem(status))

        if result.error:
            return

        self._export_btn.setEnabled(True)

        # Append top results to live table
        hits = result.single_results or result.double_results
        best_ratio = 0.0
        best_mod = ""
        for r in hits[:5]:  # Show top 5 per variant to keep table manageable
            table_row = self._results_table.rowCount()
            self._results_table.insertRow(table_row)
            if result.mode == "single":
                mod_str = r.modification
                d21, d22, ratio = r.delta_dc21, r.delta_dc22, r.delta_ratio
            else:
                mod_str = f"{r.mod_1}@{r.position_1}+{r.mod_2}@{r.position_2}"
                d21, d22, ratio = r.delta_dc21, r.delta_dc22, r.delta_ratio
                pos_str = f"{r.position_1},{r.position_2}"
            self._results_table.setItem(table_row, 0, QTableWidgetItem(result.variant_id))
            pos_display = str(getattr(r, "position", "—"))
            self._results_table.setItem(table_row, 1, QTableWidgetItem(pos_display))
            self._results_table.setItem(table_row, 2, QTableWidgetItem(mod_str))
            self._results_table.setItem(table_row, 3, QTableWidgetItem(f"{d21:+.4f}"))
            self._results_table.setItem(table_row, 4, QTableWidgetItem(f"{d22:+.4f}"))
            self._results_table.setItem(table_row, 5, QTableWidgetItem(f"{ratio:+.4f}"))
            if abs(ratio) > abs(best_ratio):
                best_ratio = ratio
                best_mod = mod_str

        # Upsert summary row
        if hits:
            self._upsert_summary(result.variant_id, best_mod, best_ratio, result.mode)

    def _on_progress(self, n_done: int, n_total: int) -> None:
        self._progress_label.setText(f"{n_done} / {n_total}")
        pct = int(n_done / n_total * 100) if n_total else 0
        self._progress_bar.setValue(pct)

    def _on_batch_done(self) -> None:
        self._run_btn.setEnabled(True)
        self._cancel_btn.setEnabled(False)
        self._progress_bar.setValue(100)

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _on_run(self) -> None:
        if not self._pending_jobs:
            return
        jobs = list(self._pending_jobs)
        self._run_btn.setEnabled(False)
        self._cancel_btn.setEnabled(True)
        self._progress_bar.setValue(0)
        self._progress_label.setText(f"0 / {len(jobs)}")
        self.run_requested.emit(jobs)

    def _add_group(self, group: str) -> None:
        if self._dataset is None:
            return
        mode = self._mode_combo.currentText()
        for vi in self._dataset.variants:
            if vi.group == group:
                self.add_job(vi.variant, mode)

    def _on_jump(self) -> None:
        # Check results_table first, then summary_table
        for table in (self._results_table, self._summary_table):
            selected = table.selectedItems()
            if selected:
                row = selected[0].row()
                vid = table.item(row, 0)
                if vid:
                    self.jump_to_variant.emit(vid.text())
                    return

    def _upsert_summary(
        self, variant_id: str, best_mod: str, best_ratio: float, mode: str
    ) -> None:
        existing_row = self._summary_row.get(variant_id)
        if existing_row is not None:
            # Only update if this result is better
            if abs(best_ratio) <= abs(self._summary_best.get(variant_id, 0.0)):
                return
            self._summary_table.setItem(existing_row, 1, QTableWidgetItem(best_mod))
            self._summary_table.setItem(
                existing_row, 2, QTableWidgetItem(f"{best_ratio:+.4f}")
            )
        else:
            row = self._summary_table.rowCount()
            self._summary_table.insertRow(row)
            self._summary_table.setItem(row, 0, QTableWidgetItem(variant_id))
            self._summary_table.setItem(row, 1, QTableWidgetItem(best_mod))
            self._summary_table.setItem(row, 2, QTableWidgetItem(f"{best_ratio:+.4f}"))
            self._summary_table.setItem(row, 3, QTableWidgetItem(mode))
            self._summary_row[variant_id] = row
        self._summary_best[variant_id] = best_ratio

    def _export_csv(self) -> None:
        """Export all results and summary to a CSV file."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_name = f"batch_results_{timestamp}.csv"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Batch Results", default_name, "CSV (*.csv)"
        )
        if not path:
            return

        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            # Detailed results
            writer.writerow(["Section: Detailed Results"])
            writer.writerow(["Variant", "Position", "Modification",
                             "Delta_DC21", "Delta_DC22", "Delta_Ratio"])
            for row in range(self._results_table.rowCount()):
                writer.writerow([
                    self._results_table.item(row, col).text()
                    if self._results_table.item(row, col) else ""
                    for col in range(6)
                ])
            writer.writerow([])
            # Summary
            writer.writerow(["Section: Best Hit per Variant"])
            writer.writerow(["Variant", "Best_Mod", "Delta_Ratio", "Mode"])
            for row in range(self._summary_table.rowCount()):
                writer.writerow([
                    self._summary_table.item(row, col).text()
                    if self._summary_table.item(row, col) else ""
                    for col in range(4)
                ])
