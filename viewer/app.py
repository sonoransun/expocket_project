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
from viewer.data.mock_data import enrich_dataset
from viewer.data.schema import EnrichedVariantDataset, VariantDataset
from viewer.interaction.controller import InteractionController
from viewer.landscape.widgets import DataLandscapeWidget
from viewer.rna3d.pocket_scene import DicerPocketOverlay
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
            "QSpinBox { background: #16213e; color: #e0e0e0; padding: 4px; }"
            "QTableWidget { background: #1e1e1e; color: #e0e0e0; "
            "gridline-color: #333; }"
            "QHeaderView::section { background: #16213e; color: #e0e0e0; "
            "padding: 4px; border: 1px solid #333; }"
        )

        # Load and enrich data
        base_dataset = dataset or load_dataset()
        self.dataset = enrich_dataset(base_dataset)

        # Central area: horizontal splitter with two 3D panels
        splitter = QSplitter(Qt.Orientation.Horizontal)

        self.rna_widget = RNAViewWidget()
        splitter.addWidget(self.rna_widget)

        self.landscape_widget = DataLandscapeWidget()
        splitter.addWidget(self.landscape_widget)

        splitter.setSizes([600, 1000])
        self.setCentralWidget(splitter)

        # DICER pocket overlay
        self._pocket_overlay = DicerPocketOverlay()
        self.rna_widget.scene_manager.scene.add(self._pocket_overlay.group)

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

        # 2D Analysis panel (bottom dock)
        self._init_2d_panel()

        # Chemistry panel (right dock)
        self._init_chemistry_panel()

        # Gene panel (bottom dock)
        self._init_gene_panel()

        # Synthesis panel (right dock)
        self._init_synthesis_panel()

        # Replacement panel (right dock)
        self._init_replacement_panel()

        # Batch screening panel (bottom dock)
        self._init_batch_panel()

        # Gene editing panel (right dock)
        self._init_editing_panel()

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

    def _init_2d_panel(self) -> None:
        """Initialize the 2D analysis tab widget."""
        try:
            from viewer.views2d.widgets import AnalysisTabWidget

            self.analysis_tabs = AnalysisTabWidget()
            dock = QDockWidget("2D Analysis", self)
            dock.setWidget(self.analysis_tabs)
            self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)
            self.analysis_tabs.set_dataset(self.dataset)
        except ImportError:
            self.analysis_tabs = None

    def _init_chemistry_panel(self) -> None:
        """Initialize the chemistry lab dock."""
        try:
            from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
            from viewer.chemistry.modification_engine import ModificationEngine
            from viewer.chemistry.virtual_screen import VirtualScreener
            from viewer.ui.chemistry_panel import ChemistryPanel

            self._mod_engine = ModificationEngine(self.dataset)
            self._predictor = CleavageSitePredictor(self.dataset)
            self._screener = VirtualScreener(
                self.dataset, self._mod_engine, self._predictor
            )

            self.chemistry_panel = ChemistryPanel()
            self.chemistry_panel.set_engines(
                self._mod_engine, self._screener, self.dataset
            )
            dock = QDockWidget("Chemistry Lab", self)
            dock.setWidget(self.chemistry_panel)
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        except ImportError:
            self.chemistry_panel = None

    def _init_gene_panel(self) -> None:
        """Initialize the gene targets dock."""
        try:
            from viewer.chemistry.mirna_context import MiRNAContextAnalyzer
            from viewer.ui.gene_panel import GenePanel

            self._mirna_analyzer = MiRNAContextAnalyzer()
            self.gene_panel = GenePanel()
            self.gene_panel.set_analyzer(self._mirna_analyzer)
            dock = QDockWidget("Gene Targets", self)
            dock.setWidget(self.gene_panel)
            self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)
        except ImportError:
            self.gene_panel = None

    def _init_synthesis_panel(self) -> None:
        """Initialize the synthesis planning dock."""
        try:
            from viewer.chemistry.synthesis_pathway import SynthesisPlanner
            from viewer.ui.synthesis_panel import SynthesisPanel

            self._synthesis_planner = SynthesisPlanner(self.dataset)
            self.synthesis_panel = SynthesisPanel()
            self.synthesis_panel.set_planner(self._synthesis_planner, self.dataset)
            dock = QDockWidget("Synthesis Plan", self)
            dock.setWidget(self.synthesis_panel)
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        except ImportError:
            self.synthesis_panel = None
            self._synthesis_planner = None

    def _init_replacement_panel(self) -> None:
        """Initialize the replacement chemistry dock."""
        try:
            from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
            from viewer.chemistry.double_screen import DoubleReplacementScreener
            from viewer.chemistry.modification_engine import ModificationEngine
            from viewer.ui.replacement_panel import ReplacementPanel

            # Reuse existing engines if chemistry panel created them
            if not hasattr(self, '_mod_engine'):
                self._mod_engine = ModificationEngine(self.dataset)
            if not hasattr(self, '_predictor'):
                self._predictor = CleavageSitePredictor(self.dataset)

            self._double_screener = DoubleReplacementScreener(
                self.dataset, self._predictor, self._mod_engine
            )
            self.replacement_panel = ReplacementPanel()
            self.replacement_panel.set_engines(
                self._mod_engine, self._predictor, self._double_screener, self.dataset
            )
            dock = QDockWidget("Replacement Chemistry", self)
            dock.setWidget(self.replacement_panel)
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        except ImportError:
            self.replacement_panel = None
            self._double_screener = None

    def _init_batch_panel(self) -> None:
        """Initialize the batch screening dock."""
        try:
            from viewer.chemistry.batch_screener import BatchScreener
            from viewer.ui.batch_panel import BatchPanel

            self._batch_screener = BatchScreener(
                screener=self._screener,
                double_screener=self._double_screener,
                dataset=self.dataset,
            )
            self.batch_panel = BatchPanel()
            self.batch_panel.set_dataset(self.dataset)
            dock = QDockWidget("Batch Screening", self)
            dock.setWidget(self.batch_panel)
            self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)
            self._batch_worker = None
        except (ImportError, AttributeError):
            self.batch_panel = None
            self._batch_worker = None

    def _connect_signals(self) -> None:
        # Sidebar -> controller
        self.sidebar.cleavage_site_changed.connect(
            self.controller.change_cleavage_site
        )
        self.sidebar.enzyme_changed.connect(self._on_enzyme_changed)
        self.sidebar.color_mode_changed.connect(
            self.controller.change_color_mode
        )
        self.sidebar.landscape_mode_changed.connect(
            self.controller.change_landscape_mode
        )
        self.sidebar.pocket_toggled.connect(self._on_pocket_toggled)

        # Controller -> info panel
        self.controller.variant_selected.connect(self._on_variant_selected)

        # Controller -> 2D analysis
        if self.analysis_tabs is not None:
            self.controller.cleavage_site_changed.connect(
                self.analysis_tabs.set_cleavage_site
            )

        # Landscape clicking -> controller
        self.landscape_widget.variant_clicked.connect(
            self.controller.select_variant
        )

        # Chemistry panel signals
        if hasattr(self, 'chemistry_panel') and self.chemistry_panel is not None:
            self.chemistry_panel.modification_applied.connect(
                self.controller.on_modification_applied
            )
            self.chemistry_panel.modification_removed.connect(
                lambda vid, pos: self.controller.on_modification_applied(vid, pos, "")
            )
            self.controller.variant_selected.connect(
                self.chemistry_panel.set_variant
            )

        # Synthesis panel signals
        if hasattr(self, 'synthesis_panel') and self.synthesis_panel is not None:
            self.controller.variant_selected.connect(
                self.synthesis_panel.set_variant
            )
            self.controller.synthesis_updated.connect(
                self._on_synthesis_updated
            )
            self.synthesis_panel.synthesis_updated.connect(
                self._on_synthesis_plan_ready
            )

        # Replacement panel signals
        if hasattr(self, 'replacement_panel') and self.replacement_panel is not None:
            self.controller.variant_selected.connect(
                self.replacement_panel.set_variant
            )
            self.replacement_panel.replacement_applied.connect(
                self.controller.on_replacement_applied
            )
            self.replacement_panel.comparison_ready.connect(
                self._on_comparison_ready
            )

        # 2D analysis: wire chemistry engines and variant updates
        if self.analysis_tabs is not None:
            if hasattr(self, '_predictor') and hasattr(self, '_mod_engine'):
                self.analysis_tabs.set_chemistry_engines(
                    self._predictor, self._mod_engine
                )
            self.controller.variant_selected.connect(
                self.analysis_tabs.set_current_variant
            )

        # Gene editing panel signals
        if hasattr(self, 'editing_panel') and self.editing_panel is not None:
            self.controller.variant_selected.connect(self.editing_panel.set_variant)
            self.editing_panel.edit_designed.connect(self._on_edit_designed)
            self.editing_panel.edit_applied.connect(self._on_edit_applied)
            if self.analysis_tabs is not None and hasattr(self.analysis_tabs, 'editing_site_clicked'):
                self.analysis_tabs.editing_site_clicked.connect(
                    self.editing_panel.select_site
                )

        # Batch screening panel
        if hasattr(self, 'batch_panel') and self.batch_panel is not None:
            self.batch_panel.run_requested.connect(self._on_batch_run)
            self.batch_panel.cancel_requested.connect(self._on_batch_cancel)
            self.batch_panel.jump_to_variant.connect(self.controller.select_variant)
            if hasattr(self, 'chemistry_panel') and self.chemistry_panel is not None:
                self.chemistry_panel.add_to_batch_requested.connect(
                    self.batch_panel.add_job
                )
            if hasattr(self, 'replacement_panel') and self.replacement_panel is not None:
                self.replacement_panel.add_to_batch_requested.connect(
                    self.batch_panel.add_job
                )

    def _on_variant_selected(self, variant_id: str) -> None:
        self.info_panel.update_variant(variant_id, self.controller.dataset)

        # Update gene panel with cleavage data
        if hasattr(self, 'gene_panel') and self.gene_panel is not None:
            rec21 = self.dataset.get_cleavage(variant_id, 21)
            rec22 = self.dataset.get_cleavage(variant_id, 22)
            dc21 = rec21.mean_accuracy if rec21 else 0.0
            dc22 = rec22.mean_accuracy if rec22 else 0.0
            self.gene_panel.update_cleavage(dc21, dc22)

        # Build pocket overlay if layout is available
        if self.dataset.dicer_pocket and self.rna_widget.scene_manager._layout:
            self._pocket_overlay.build(
                self.dataset.dicer_pocket,
                self.rna_widget.scene_manager._layout,
            )

        # Auto-populate editing map for all tools
        if (
            hasattr(self, "_target_finder")
            and self._target_finder is not None
            and self.analysis_tabs is not None
        ):
            try:
                sites = self._target_finder.find_all_tools(variant_id, self.dataset)
                self._last_editing_sites = sites
                self.analysis_tabs.update_editing_sites(variant_id, self.dataset, sites)
            except Exception:
                pass

    def _on_enzyme_changed(self, enzyme: str) -> None:
        # Cancel any in-flight batch before replacing the dataset/engines
        if hasattr(self, '_batch_worker') and self._batch_worker and self._batch_worker.isRunning():
            self._batch_worker.cancel()
            self._batch_worker.wait()

        base_dataset = load_dataset(enzyme=enzyme)
        self.dataset = enrich_dataset(base_dataset)
        self.controller.set_dataset(self.dataset)

        if self.analysis_tabs is not None:
            self.analysis_tabs.set_dataset(self.dataset)

        # Rebuild batch screener with fresh dataset and engines
        if hasattr(self, 'batch_panel') and self.batch_panel is not None:
            try:
                from viewer.chemistry.batch_screener import BatchScreener
                self._batch_screener = BatchScreener(
                    screener=self._screener,
                    double_screener=self._double_screener,
                    dataset=self.dataset,
                )
                self.batch_panel.set_dataset(self.dataset)
            except (ImportError, AttributeError):
                pass

    def _on_synthesis_updated(self, variant_id: str) -> None:
        """Refresh synthesis panel when modifications change."""
        if hasattr(self, 'synthesis_panel') and self.synthesis_panel is not None:
            self.synthesis_panel.refresh()

    def _on_synthesis_plan_ready(self, variant_id: str) -> None:
        """Push synthesis plan to 2D diagram tab."""
        if self.analysis_tabs is not None and hasattr(self, 'synthesis_panel') and self.synthesis_panel is not None:
            self.analysis_tabs.update_synthesis(self.synthesis_panel.current_plan)

    def _on_comparison_ready(self, comparison: dict) -> None:
        """Push replacement comparison data to 2D tab."""
        if self.analysis_tabs is not None:
            self.analysis_tabs.update_replacement(comparison)

    def _on_batch_run(self, jobs: list) -> None:
        from viewer.workers.batch_worker import BatchWorker

        if self._batch_worker and self._batch_worker.isRunning():
            return
        self._batch_worker = BatchWorker(self._batch_screener, jobs)
        self._batch_worker.result_ready.connect(self.batch_panel._on_result)
        self._batch_worker.progress_updated.connect(self.batch_panel._on_progress)
        self._batch_worker.batch_finished.connect(self.batch_panel._on_batch_done)
        self._batch_worker.start()

    def _on_batch_cancel(self) -> None:
        if self._batch_worker:
            self._batch_worker.cancel()

    def _on_pocket_toggled(self, visible: bool) -> None:
        self._pocket_overlay.set_visible(visible)

    def _init_editing_panel(self) -> None:
        """Initialize the gene editing dock panel."""
        try:
            from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
            from viewer.genome_editing.edit_engine import EditEngine
            from viewer.genome_editing.guide_designer import GuideDesigner
            from viewer.genome_editing.impact_predictor import EditImpactPredictor
            from viewer.genome_editing.off_target import OffTargetScorer
            from viewer.genome_editing.target_finder import TargetFinder
            from viewer.ui.editing_panel import EditingPanel

            if not hasattr(self, "_predictor"):
                self._predictor = CleavageSitePredictor(self.dataset)

            self._target_finder = TargetFinder()
            self._guide_designer = GuideDesigner()
            self._edit_engine = EditEngine()
            self._edit_impact = EditImpactPredictor(self._predictor, self.dataset)
            self._off_target_scorer = OffTargetScorer()
            self._last_editing_sites: dict = {}

            self.editing_panel = EditingPanel(
                dataset=self.dataset,
                target_finder=self._target_finder,
                guide_designer=self._guide_designer,
                edit_engine=self._edit_engine,
                impact_predictor=self._edit_impact,
                off_target_scorer=self._off_target_scorer,
                parent=self,
            )
            dock = QDockWidget("Gene Editing", self)
            dock.setWidget(self.editing_panel)
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        except Exception as exc:
            print(f"[editing] panel not available: {exc}")
            self.editing_panel = None

    def _on_edit_designed(self, variant_id: str, design: object) -> None:
        """Push editing sites to the 2D map when a design is ready."""
        from viewer.genome_editing.target_finder import TargetFinder
        if not hasattr(self, "_target_finder"):
            return
        try:
            sites = self._target_finder.find_all_tools(variant_id, self.dataset)
            self._last_editing_sites = sites
            if self.analysis_tabs is not None:
                self.analysis_tabs.update_editing_sites(variant_id, self.dataset, sites)
        except Exception:
            pass

    def _on_edit_applied(self, variant_id: str, modified_rna: str) -> None:
        """Preview a modified RNA sequence in the 3D RNA viewer."""
        from viewer.data.schema import VariantInfo
        variant = self.dataset.get_variant(variant_id)
        if variant is None:
            return
        # Construct a temporary VariantInfo with the modified sequence
        temp = VariantInfo(
            variant=f"{variant_id}_edited",
            group=variant.group,
            randomized_nts=variant.randomized_nts,
            pre_mirna_sequence=modified_rna,
            concrete_struct=variant.concrete_struct,
            flanking_length_5p=variant.flanking_length_5p,
        )
        cleavage_data = self.dataset.cleavage_data.get(variant_id, [])
        self.rna_widget.scene_manager.set_variant(temp, cleavage_data, None)
