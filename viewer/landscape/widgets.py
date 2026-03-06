"""Qt widget wrapping the pygfx canvas for the data landscape panel."""

from __future__ import annotations

import PySide6  # noqa: F401 — must import before rendercanvas.qt
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QVBoxLayout, QWidget
from rendercanvas.qt import QRenderWidget
import pygfx

from viewer.data.schema import VariantDataset
from viewer.landscape.scene import DataLandscapeScene


class DataLandscapeWidget(QWidget):
    """Qt widget displaying the 3D data landscape of 256 variants."""

    variant_clicked = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.canvas = QRenderWidget(parent=self)
        layout.addWidget(self.canvas)

        self.renderer = pygfx.renderers.WgpuRenderer(self.canvas)
        self.scene_manager = DataLandscapeScene()

        self.controller = pygfx.OrbitController()
        self.controller.add_camera(self.scene_manager.camera)
        self.controller.register_events(self.renderer)

        self.canvas.request_draw(self._animate)

    def _animate(self) -> None:
        self.renderer.render(
            self.scene_manager.scene, self.scene_manager.camera
        )
        self.canvas.request_draw(self._animate)

    def load_dataset(
        self, dataset: VariantDataset, cleavage_site: int = 21
    ) -> None:
        """Build the scatter plot for the given dataset."""
        self.scene_manager.build_scatter(dataset, cleavage_site)
