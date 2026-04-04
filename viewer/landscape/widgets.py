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

        # Wire picking via double-click to avoid interfering with orbit
        self.renderer.add_event_handler(self._on_pointer_down, "double_click")

    def _animate(self) -> None:
        self.renderer.render(
            self.scene_manager.scene, self.scene_manager.camera
        )
        self.canvas.request_draw(self._animate)

    def _on_pointer_down(self, event) -> None:
        """Handle double-click for variant picking."""
        try:
            # Use pygfx picking if available
            pick_info = self.renderer.get_pick_info(event)
            if pick_info and "world_object" in pick_info:
                mesh = pick_info["world_object"]
                if isinstance(mesh, pygfx.Mesh):
                    vid = self.scene_manager.get_variant_for_mesh(mesh)
                    if vid:
                        self.variant_clicked.emit(vid)
                        return
        except (AttributeError, TypeError) as exc:
            import logging
            logging.getLogger(__name__).debug("Picking via pygfx failed: %s", exc)

        # Fallback: use nearest-variant by world position
        try:
            if hasattr(event, "x") and hasattr(event, "y"):
                import numpy as np
                # Approximate: project screen center to world
                pos = np.array([event.x, event.y, 0.0])
                vid = self.scene_manager.get_nearest_variant(pos)
                if vid:
                    self.variant_clicked.emit(vid)
        except Exception as exc:
            import logging
            logging.getLogger(__name__).debug("Fallback picking failed: %s", exc)

    def load_dataset(
        self, dataset: VariantDataset, cleavage_site: int = 21
    ) -> None:
        """Build the scatter plot for the given dataset."""
        self.scene_manager.build_scatter(dataset, cleavage_site)
