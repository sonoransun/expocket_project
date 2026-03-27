"""QTabWidget hosting all 2D analysis views with interactive event handling."""

from __future__ import annotations

from typing import TYPE_CHECKING

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QTabWidget, QVBoxLayout, QWidget

from viewer.views2d.contact_map import plot_contact_map
from viewer.views2d.modification_impact import clear_impact_cache, plot_modification_impact
from viewer.views2d.property_heatmap import plot_property_heatmap
from viewer.views2d.sar_matrix import plot_sar_matrix

if TYPE_CHECKING:
    from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
    from viewer.chemistry.modification_engine import ModificationEngine
    from viewer.data.schema import EnrichedVariantDataset, SynthesisPlan
    from viewer.views2d.editing_map import EditingMapWidget


class _CanvasTab(QWidget):
    """A single tab holding a matplotlib FigureCanvas with event support."""

    cell_clicked = Signal(int, int)
    cell_hovered = Signal(int, int, str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.figure.set_facecolor("#1e1e1e")
        self.figure.patch.set_facecolor("#1e1e1e")
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)
        self._annotation = None
        self._data_matrix = None  # set externally for tooltip lookup
        self._row_labels: list[str] = []
        self._col_labels: list[str] = []

    def connect_events(self) -> None:
        """Wire matplotlib event handlers for hover and click."""
        self.canvas.mpl_connect("motion_notify_event", self._on_hover)
        self.canvas.mpl_connect("button_press_event", self._on_click)

    def set_data_for_tooltips(
        self,
        matrix,
        row_labels: list[str] | None = None,
        col_labels: list[str] | None = None,
    ) -> None:
        """Store data matrix for tooltip display."""
        self._data_matrix = matrix
        self._row_labels = row_labels or []
        self._col_labels = col_labels or []

    def _on_hover(self, event) -> None:
        if event.inaxes is None or self._data_matrix is None:
            if self._annotation is not None:
                self._annotation.set_visible(False)
                self.canvas.draw_idle()
            return
        ax = event.inaxes
        col = int(round(event.xdata))
        row = int(round(event.ydata))
        import numpy as np
        data = np.asarray(self._data_matrix)
        if 0 <= row < data.shape[0] and 0 <= col < data.shape[1]:
            val = data[row, col]
            rlabel = self._row_labels[row] if row < len(self._row_labels) else str(row)
            clabel = self._col_labels[col] if col < len(self._col_labels) else str(col)
            text = f"{rlabel} / {clabel}\n{val:.4f}"
            if self._annotation is None:
                self._annotation = ax.annotate(
                    text, xy=(col, row), fontsize=6,
                    bbox=dict(boxstyle="round,pad=0.3", fc="#333", ec="#666", alpha=0.9),
                    color="#eee", ha="left", va="bottom",
                    xytext=(8, 8), textcoords="offset points",
                )
            else:
                self._annotation.set_text(text)
                self._annotation.xy = (col, row)
                self._annotation.set_visible(True)
            self.canvas.draw_idle()
        else:
            if self._annotation is not None:
                self._annotation.set_visible(False)
                self.canvas.draw_idle()

    def _on_click(self, event) -> None:
        if event.inaxes is None:
            return
        col = int(round(event.xdata))
        row = int(round(event.ydata))
        self.cell_clicked.emit(row, col)

    def redraw(self) -> None:
        self.canvas.draw_idle()


class AnalysisTabWidget(QTabWidget):
    """Tab widget containing all 2D analysis panels."""

    modification_cell_clicked = Signal(int, str)  # position, mod_code
    editing_site_clicked = Signal(object)          # emits TargetSite

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        self._heatmap_tab = _CanvasTab()
        self._sar_tab = _CanvasTab()
        self._contact_tab = _CanvasTab()
        self._mod_impact_tab = _CanvasTab()
        self._synthesis_tab = _CanvasTab()
        self._replacement_tab = _CanvasTab()

        self.addTab(self._heatmap_tab, "Properties")
        self.addTab(self._sar_tab, "SAR")
        self.addTab(self._contact_tab, "Contacts")
        self.addTab(self._mod_impact_tab, "Mod Impact")
        self.addTab(self._synthesis_tab, "Synthesis")
        self.addTab(self._replacement_tab, "Replacements")

        # Gene editing map tab
        self._editing_map: "EditingMapWidget | None" = None
        try:
            from viewer.views2d.editing_map import EditingMapWidget
            self._editing_map = EditingMapWidget()
            self.addTab(self._editing_map, "Editing")
            self._editing_map.site_clicked.connect(self.editing_site_clicked)
        except Exception:
            self._editing_map = None

        self._dataset: EnrichedVariantDataset | None = None
        self._cleavage_site: int = 21
        self._predictor: CleavageSitePredictor | None = None
        self._engine: ModificationEngine | None = None
        self._current_variant: str | None = None

        # Wire click events
        self._mod_impact_tab.connect_events()
        self._mod_impact_tab.cell_clicked.connect(self._on_mod_impact_click)
        self._heatmap_tab.connect_events()
        self._sar_tab.connect_events()
        self._contact_tab.connect_events()

        # Style dark theme
        self.setStyleSheet("""
            QTabWidget::pane { border: 1px solid #444; background: #1e1e1e; }
            QTabBar::tab { background: #2d2d2d; color: #ccc; padding: 6px 14px;
                           border: 1px solid #444; border-bottom: none; }
            QTabBar::tab:selected { background: #1e1e1e; color: #fff; }
        """)

    def set_chemistry_engines(
        self,
        predictor: CleavageSitePredictor,
        engine: ModificationEngine,
    ) -> None:
        """Set chemistry engines for real modification impact computation."""
        self._predictor = predictor
        self._engine = engine

    def set_current_variant(self, variant_id: str) -> None:
        """Update the current variant and refresh variant-specific views."""
        self._current_variant = variant_id
        self._refresh_mod_impact()

    def set_dataset(self, dataset: EnrichedVariantDataset) -> None:
        self._dataset = dataset
        clear_impact_cache()
        self.refresh_all()

    def set_cleavage_site(self, site: int) -> None:
        self._cleavage_site = site
        self._refresh_sar()

    def refresh_all(self) -> None:
        self._refresh_heatmap()
        self._refresh_sar()
        self._refresh_contacts()
        self._refresh_mod_impact()

    def update_synthesis(self, plan: SynthesisPlan | None) -> None:
        """Update the synthesis diagram tab."""
        fig = self._synthesis_tab.figure
        try:
            from viewer.views2d.synthesis_diagram import plot_synthesis_diagram
            plot_synthesis_diagram(fig, plan)
        except Exception:
            fig.clear()
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, "No synthesis plan", ha="center", va="center", color="#ccc")
        self._apply_dark_axes(fig)
        self._synthesis_tab.redraw()

    def update_editing_sites(
        self,
        variant_id: str,
        dataset: "EnrichedVariantDataset",
        sites_by_tool: dict,
    ) -> None:
        """Refresh the editing map tab with new site data."""
        if self._editing_map is not None:
            self._editing_map.update_variant(variant_id, dataset, sites_by_tool)

    def update_replacement(self, comparison: dict | None) -> None:
        """Update the replacement comparison tab."""
        fig = self._replacement_tab.figure
        try:
            from viewer.views2d.replacement_comparison import plot_replacement_comparison
            if comparison is not None:
                plot_replacement_comparison(
                    fig,
                    comparison.get("original", {}),
                    comparison.get("single_1"),
                    comparison.get("single_2"),
                    comparison.get("double"),
                )
            else:
                fig.clear()
                ax = fig.add_subplot(111)
                ax.text(0.5, 0.5, "Select replacements to compare", ha="center", va="center", color="#ccc")
        except Exception:
            fig.clear()
        self._apply_dark_axes(fig)
        self._replacement_tab.redraw()

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
        plot_modification_impact(
            fig,
            seq_length=seq_len,
            variant_id=self._current_variant,
            predictor=self._predictor,
            engine=self._engine,
        )
        self._apply_dark_axes(fig)
        self._mod_impact_tab.redraw()

        # Set tooltip data for the mod impact grid
        from viewer.encoding.modification_db import MODIFICATIONS_DB
        from viewer.views2d.modification_impact import _impact_cache
        mod_codes = list(MODIFICATIONS_DB.keys())
        if self._current_variant and self._current_variant in _impact_cache:
            import numpy as np
            matrix = _impact_cache[self._current_variant]
            # The image is transposed (mods on Y, positions on X)
            self._mod_impact_tab.set_data_for_tooltips(
                matrix.T, row_labels=mod_codes,
                col_labels=[str(i) for i in range(seq_len)],
            )

    def _on_mod_impact_click(self, row: int, col: int) -> None:
        """Handle click on modification impact grid cell."""
        from viewer.encoding.modification_db import MODIFICATIONS_DB
        mod_codes = list(MODIFICATIONS_DB.keys())
        if 0 <= row < len(mod_codes):
            self.modification_cell_clicked.emit(col, mod_codes[row])

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
