"""QTabWidget hosting all four 2D analysis views."""

from __future__ import annotations

from typing import TYPE_CHECKING

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtWidgets import QTabWidget, QVBoxLayout, QWidget

from viewer.views2d.contact_map import plot_contact_map
from viewer.views2d.modification_impact import plot_modification_impact
from viewer.views2d.property_heatmap import plot_property_heatmap
from viewer.views2d.sar_matrix import plot_sar_matrix

if TYPE_CHECKING:
    from viewer.data.schema import EnrichedVariantDataset
    from viewer.encoding.protein_descriptors import DicerPocketModel


class _CanvasTab(QWidget):
    """A single tab holding a matplotlib FigureCanvas."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.figure.set_facecolor("#1e1e1e")
        self.figure.patch.set_facecolor("#1e1e1e")
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

    def redraw(self) -> None:
        self.canvas.draw_idle()


class AnalysisTabWidget(QTabWidget):
    """Tab widget containing all 2D analysis panels."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self._heatmap_tab = _CanvasTab()
        self._sar_tab = _CanvasTab()
        self._contact_tab = _CanvasTab()
        self._mod_impact_tab = _CanvasTab()

        self.addTab(self._heatmap_tab, "Properties")
        self.addTab(self._sar_tab, "SAR")
        self.addTab(self._contact_tab, "Contacts")
        self.addTab(self._mod_impact_tab, "Mod Impact")

        self._dataset: EnrichedVariantDataset | None = None
        self._cleavage_site: int = 21

        # Style dark theme
        self.setStyleSheet("""
            QTabWidget::pane { border: 1px solid #444; background: #1e1e1e; }
            QTabBar::tab { background: #2d2d2d; color: #ccc; padding: 6px 14px;
                           border: 1px solid #444; border-bottom: none; }
            QTabBar::tab:selected { background: #1e1e1e; color: #fff; }
        """)

    def set_dataset(self, dataset: EnrichedVariantDataset) -> None:
        self._dataset = dataset
        self.refresh_all()

    def set_cleavage_site(self, site: int) -> None:
        self._cleavage_site = site
        self._refresh_sar()

    def refresh_all(self) -> None:
        self._refresh_heatmap()
        self._refresh_sar()
        self._refresh_contacts()
        self._refresh_mod_impact()

    def _refresh_heatmap(self) -> None:
        if self._dataset is None:
            return
        fig = self._heatmap_tab.figure
        plot_property_heatmap(fig, self._dataset)
        self._apply_dark_axes(fig)
        self._heatmap_tab.redraw()

    def _refresh_sar(self) -> None:
        if self._dataset is None:
            return
        fig = self._sar_tab.figure
        plot_sar_matrix(fig, self._dataset, self._cleavage_site)
        self._apply_dark_axes(fig)
        self._sar_tab.redraw()

    def _refresh_contacts(self) -> None:
        if self._dataset is None:
            return
        pocket = self._dataset.dicer_pocket
        seq_len = 63
        if self._dataset.variants:
            seq_len = len(self._dataset.variants[0].pre_mirna_sequence)
        fig = self._contact_tab.figure
        plot_contact_map(fig, pocket, seq_len)
        self._apply_dark_axes(fig)
        self._contact_tab.redraw()

    def _refresh_mod_impact(self) -> None:
        fig = self._mod_impact_tab.figure
        seq_len = 63
        if self._dataset and self._dataset.variants:
            seq_len = len(self._dataset.variants[0].pre_mirna_sequence)
        plot_modification_impact(fig, seq_length=seq_len)
        self._apply_dark_axes(fig)
        self._mod_impact_tab.redraw()

    @staticmethod
    def _apply_dark_axes(fig: Figure) -> None:
        """Apply dark theme to all axes in the figure."""
        for ax in fig.get_axes():
            ax.set_facecolor("#2d2d2d")
            ax.tick_params(colors="#ccc", labelsize=6)
            ax.xaxis.label.set_color("#ccc")
            ax.yaxis.label.set_color("#ccc")
            ax.title.set_color("#ccc")
            for spine in ax.spines.values():
                spine.set_color("#555")
        if fig._suptitle is not None:
            fig._suptitle.set_color("#eee")
