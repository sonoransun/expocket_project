"""Sidebar control panel with filters, selectors, and toggles."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QGroupBox,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from viewer.config import GROUPS


class SidebarWidget(QWidget):
    """Control panel for filtering and display options."""

    cleavage_site_changed = Signal(int)
    enzyme_changed = Signal(str)
    group_filter_changed = Signal(list)
    color_mode_changed = Signal(str)
    landscape_mode_changed = Signal(str)
    pocket_toggled = Signal(bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setSpacing(12)

        # Cleavage site selector
        layout.addWidget(QLabel("Cleavage Site:"))
        self.site_selector = QComboBox()
        self.site_selector.addItems(["DC21", "DC22", "DC20", "DC23"])
        self.site_selector.currentTextChanged.connect(self._on_site_changed)
        layout.addWidget(self.site_selector)

        # Enzyme selector
        layout.addWidget(QLabel("Enzyme:"))
        self.enzyme_selector = QComboBox()
        self.enzyme_selector.addItems(["Human DICER", "Fly DCR-1"])
        self.enzyme_selector.currentTextChanged.connect(self._on_enzyme_changed)
        layout.addWidget(self.enzyme_selector)

        # Group filter
        group_box = QGroupBox("5' Nucleotide Group")
        group_layout = QVBoxLayout(group_box)
        self._group_checks: dict[str, QCheckBox] = {}
        for g in GROUPS:
            cb = QCheckBox(g)
            cb.setChecked(True)
            cb.stateChanged.connect(self._on_group_filter)
            group_layout.addWidget(cb)
            self._group_checks[g] = cb
        layout.addWidget(group_box)

        # Annotation toggles
        ann_box = QGroupBox("Annotations")
        ann_layout = QVBoxLayout(ann_box)
        self.show_cleavage = QCheckBox("Cleavage planes")
        self.show_cleavage.setChecked(True)
        ann_layout.addWidget(self.show_cleavage)
        self.show_randomized = QCheckBox("Highlight randomized nts")
        self.show_randomized.setChecked(True)
        ann_layout.addWidget(self.show_randomized)
        self.show_pairs = QCheckBox("Base pair lines")
        self.show_pairs.setChecked(True)
        ann_layout.addWidget(self.show_pairs)
        layout.addWidget(ann_box)

        # Encoding section
        enc_box = QGroupBox("Encoding / Coloring")
        enc_layout = QVBoxLayout(enc_box)
        enc_layout.addWidget(QLabel("Base Color Mode:"))
        self.color_mode_combo = QComboBox()
        self.color_mode_combo.addItems([
            "nucleotide", "molecular_weight", "vdw_volume",
            "h_bond_donors", "h_bond_acceptors", "hydrophobicity_index",
            "stacking_energy_5p", "stacking_energy_3p", "sugar_pucker_preference",
        ])
        self.color_mode_combo.currentTextChanged.connect(
            lambda t: self.color_mode_changed.emit(t)
        )
        enc_layout.addWidget(self.color_mode_combo)

        self.pocket_check = QCheckBox("Show DICER pocket")
        self.pocket_check.stateChanged.connect(
            lambda s: self.pocket_toggled.emit(s == 2)
        )
        enc_layout.addWidget(self.pocket_check)
        layout.addWidget(enc_box)

        # Landscape mode
        lm_box = QGroupBox("Landscape Layout")
        lm_layout = QVBoxLayout(lm_box)
        self.landscape_mode_combo = QComboBox()
        self.landscape_mode_combo.addItems(["Grid", "PCA Property Space"])
        self.landscape_mode_combo.currentTextChanged.connect(self._on_landscape_mode)
        lm_layout.addWidget(self.landscape_mode_combo)
        layout.addWidget(lm_box)

        layout.addStretch()

    def _on_site_changed(self, text: str) -> None:
        site_map = {"DC20": 20, "DC21": 21, "DC22": 22, "DC23": 23}
        site = site_map.get(text, 21)
        self.cleavage_site_changed.emit(site)

    def _on_enzyme_changed(self, text: str) -> None:
        enzyme_map = {"Human DICER": "hdicer", "Fly DCR-1": "dcr"}
        enzyme = enzyme_map.get(text, "hdicer")
        self.enzyme_changed.emit(enzyme)

    def _on_group_filter(self) -> None:
        active = [g for g, cb in self._group_checks.items() if cb.isChecked()]
        self.group_filter_changed.emit(active)

    def _on_landscape_mode(self, text: str) -> None:
        mode_map = {"Grid": "grid", "PCA Property Space": "pca"}
        self.landscape_mode_changed.emit(mode_map.get(text, "grid"))
