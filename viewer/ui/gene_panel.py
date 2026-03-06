"""Gene targets panel showing miRNA products and downstream targets."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtWidgets import (
    QGroupBox,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from viewer.chemistry.mirna_context import MiRNAContextAnalyzer


class GenePanel(QWidget):
    """Displays miRNA product info and predicted target genes."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # Product info
        prod_box = QGroupBox("Dominant miRNA Product")
        prod_layout = QVBoxLayout(prod_box)
        self._product_label = QLabel("No data")
        self._product_label.setWordWrap(True)
        self._product_label.setStyleSheet("font-family: 'Courier New', Courier, monospace; font-size: 11px;")
        prod_layout.addWidget(self._product_label)
        self._bias_label = QLabel("")
        prod_layout.addWidget(self._bias_label)
        layout.addWidget(prod_box)

        # Target genes table
        gene_box = QGroupBox("Target Genes")
        gene_layout = QVBoxLayout(gene_box)
        self._gene_table = QTableWidget(0, 4)
        self._gene_table.setHorizontalHeaderLabels(
            ["Gene", "Binding", "Score", "Pathway"]
        )
        self._gene_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        gene_layout.addWidget(self._gene_table)
        layout.addWidget(gene_box)

        layout.addStretch()

        self._analyzer: MiRNAContextAnalyzer | None = None

    def set_analyzer(self, analyzer: MiRNAContextAnalyzer) -> None:
        self._analyzer = analyzer
        self._refresh()

    def update_cleavage(self, dc21_acc: float, dc22_acc: float) -> None:
        """Update with new cleavage accuracy values."""
        if self._analyzer is None:
            return
        # Determine dominant site
        if dc21_acc >= dc22_acc:
            self._analyzer.set_dominant_site(21)
        else:
            self._analyzer.set_dominant_site(22)

        bias = self._analyzer.get_strand_bias(dc21_acc, dc22_acc)
        self._bias_label.setText(f"<b>Strand bias:</b> {bias}")
        self._refresh()

    def _refresh(self) -> None:
        if self._analyzer is None:
            return

        product = self._analyzer.get_dominant_product()
        if product:
            self._product_label.setText(
                f"<b>{product.name}</b> ({product.strand})<br>"
                f"<code>{product.sequence}</code><br>"
                f"From DC{product.cleavage_site} cleavage"
            )
        else:
            self._product_label.setText("No dominant product")

        genes = self._analyzer.get_affected_genes()
        self._gene_table.setRowCount(0)
        for gene in genes:
            row = self._gene_table.rowCount()
            self._gene_table.insertRow(row)
            self._gene_table.setItem(row, 0, QTableWidgetItem(gene.symbol))
            self._gene_table.setItem(row, 1, QTableWidgetItem(gene.binding_type))
            self._gene_table.setItem(row, 2, QTableWidgetItem(f"{gene.context_score:.2f}"))
            self._gene_table.setItem(row, 3, QTableWidgetItem(gene.pathway))
