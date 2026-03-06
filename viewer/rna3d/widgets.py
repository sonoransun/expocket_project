"""Qt widget wrapping the pygfx canvas for the RNA hairpin panel."""

from __future__ import annotations

import PySide6  # noqa: F401 — must import before rendercanvas.qt
from PySide6.QtWidgets import QVBoxLayout, QWidget
from rendercanvas.qt import QRenderWidget
import pygfx

from viewer.data.schema import CleavageRecord, VariantInfo
from viewer.rna3d.scene import RNAHairpinScene


class RNAViewWidget(QWidget):
    """Qt widget displaying the 3D RNA hairpin structure."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.canvas = QRenderWidget(parent=self)
        layout.addWidget(self.canvas)

        self.renderer = pygfx.renderers.WgpuRenderer(self.canvas)
        self.scene_manager = RNAHairpinScene()

        self.controller = pygfx.OrbitController()
        self.controller.add_camera(self.scene_manager.camera)
        self.controller.register_events(self.renderer)

        self.canvas.request_draw(self._animate)

    def _animate(self) -> None:
        self.renderer.render(
            self.scene_manager.scene, self.scene_manager.camera
        )
        self.canvas.request_draw(self._animate)

    def set_variant(
        self,
        variant: VariantInfo,
        cleavage_data: list[CleavageRecord] | None = None,
    ) -> None:
        """Update the displayed variant."""
        self.scene_manager.set_variant(variant, cleavage_data)
