"""Info panel showing selected variant details."""

from __future__ import annotations

from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

from viewer.data.schema import CleavageRecord, VariantDataset


class InfoPanelWidget(QWidget):
    """Displays detailed information about the currently selected variant."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self._label = QLabel("Select a variant to view details.")
        self._label.setWordWrap(True)
        self._label.setStyleSheet("font-family: 'Courier New', Courier, monospace; font-size: 12px;")
        layout.addWidget(self._label)

    def update_variant(self, variant_id: str, dataset: VariantDataset) -> None:
        """Update the display for the given variant."""
        vi = dataset.get_variant(variant_id)
        if vi is None:
            self._label.setText(f"Variant {variant_id} not found.")
            return

        lines = [
            f"<b>Variant:</b> {vi.variant}",
            f"<b>Group (5' nt):</b> {vi.group}",
            f"<b>Randomized 3' nts:</b> {vi.randomized_nts}",
            f"<b>Sequence:</b> <code>{vi.pre_mirna_sequence}</code>",
            f"<b>Structure:</b> <code>{vi.concrete_struct}</code>",
            f"<b>5' flanking:</b> {vi.flanking_length_5p}",
            "",
            "<b>Cleavage Data:</b>",
        ]

        records = dataset.cleavage_data.get(variant_id, [])
        for rec in sorted(records, key=lambda r: r.cleavage_site):
            lines.append(
                f"  DC{rec.cleavage_site}: accuracy={rec.mean_accuracy:.3f}"
                f"  eff={rec.mean_positional_efficiency:.2f}"
            )

        self._label.setText("<br>".join(lines))
