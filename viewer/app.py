"""Main application window with dual-panel 3D viewer."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QMainWindow,
    QSplitter,
)

from viewer.data.loader import load_dataset
from viewer.data.schema import VariantDataset
from viewer.interaction.controller import InteractionController
from viewer.landscape.widgets import DataLandscapeWidget
from viewer.rna3d.widgets import RNAViewWidget
from viewer.ui.info_panel import InfoPanelWidget
from viewer.ui.sidebar import SidebarWidget


class MainWindow(QMainWindow):
    """Main application window with RNA 3D view and data landscape."""

    def __init__(self, dataset: VariantDataset | None = None) -> None:
        super().__init__()
        self.setWindowTitle("Pre-miRNA 3D Viewer — DICER Cleavage Analysis")
        self.resize(1600, 900)
        self.setStyleSheet(
            "QMainWindow { background-color: #1a1a2e; }"
            "QDockWidget { color: #e0e0e0; }"
            "QDockWidget::title { background: #16213e; padding: 4px; }"
            "QLabel { color: #e0e0e0; }"
            "QGroupBox { color: #e0e0e0; border: 1px solid #333; "
            "border-radius: 4px; margin-top: 8px; padding-top: 12px; }"
            "QGroupBox::title { subcontrol-origin: margin; left: 8px; }"
            "QComboBox { background: #16213e; color: #e0e0e0; padding: 4px; }"
            "QCheckBox { color: #e0e0e0; }"
            "QPushButton { background: #0f3460; color: #e0e0e0; "
            "padding: 6px 12px; border-radius: 4px; }"
        )

        # Load data
        self.dataset = dataset or load_dataset()

        # Central area: horizontal splitter with two 3D panels
        splitter = QSplitter(Qt.Orientation.Horizontal)

        self.rna_widget = RNAViewWidget()
        splitter.addWidget(self.rna_widget)

        self.landscape_widget = DataLandscapeWidget()
        splitter.addWidget(self.landscape_widget)

        splitter.setSizes([600, 1000])
        self.setCentralWidget(splitter)

        # Sidebar (right dock)
        self.sidebar = SidebarWidget()
        sidebar_dock = QDockWidget("Controls", self)
        sidebar_dock.setWidget(self.sidebar)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, sidebar_dock)

        # Info panel (bottom dock)
        self.info_panel = InfoPanelWidget()
        info_dock = QDockWidget("Variant Details", self)
        info_dock.setWidget(self.info_panel)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, info_dock)

        # Interaction controller
        self.controller = InteractionController(
            rna_scene=self.rna_widget.scene_manager,
            landscape_scene=self.landscape_widget.scene_manager,
            dataset=self.dataset,
        )

        # Wire signals
        self._connect_signals()

        # Initial state: load landscape and select first variant
        self.landscape_widget.load_dataset(self.dataset, cleavage_site=21)
        if self.dataset.variants:
            self.controller.select_variant(self.dataset.variants[0].variant)

    def _connect_signals(self) -> None:
        # Sidebar -> controller
        self.sidebar.cleavage_site_changed.connect(
            self.controller.change_cleavage_site
        )
        self.sidebar.enzyme_changed.connect(self._on_enzyme_changed)

        # Controller -> info panel
        self.controller.variant_selected.connect(self._on_variant_selected)

        # Landscape clicking -> controller
        self.landscape_widget.variant_clicked.connect(
            self.controller.select_variant
        )

    def _on_variant_selected(self, variant_id: str) -> None:
        self.info_panel.update_variant(variant_id, self.controller.dataset)

    def _on_enzyme_changed(self, enzyme: str) -> None:
        new_dataset = load_dataset(enzyme=enzyme)
        self.dataset = new_dataset
        self.controller.set_dataset(new_dataset)
