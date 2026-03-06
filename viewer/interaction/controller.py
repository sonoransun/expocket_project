"""Cross-panel interaction controller linking RNA view and data landscape."""

from __future__ import annotations

from PySide6.QtCore import QObject, Signal

from viewer.data.schema import EnrichedVariantDataset, VariantDataset
from viewer.landscape.scene import DataLandscapeScene
from viewer.rna3d.scene import RNAHairpinScene


class InteractionController(QObject):
    """Coordinates selection state between the two 3D panels."""

    variant_selected = Signal(str)
    cleavage_site_changed = Signal(int)
    enzyme_changed = Signal(str)
    modification_applied = Signal(str, int, str)  # variant_id, position, mod_code
    replacement_applied = Signal(str, dict)  # variant_id, {pos: mod_code}
    synthesis_updated = Signal(str)  # variant_id
    color_mode_changed = Signal(str)
    landscape_mode_changed = Signal(str)

    def __init__(
        self,
        rna_scene: RNAHairpinScene,
        landscape_scene: DataLandscapeScene,
        dataset: VariantDataset,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._rna = rna_scene
        self._landscape = landscape_scene
        self._dataset = dataset
        self._current_variant: str | None = None
        self._current_site: int = 21

    @property
    def current_variant(self) -> str | None:
        return self._current_variant

    @property
    def dataset(self) -> VariantDataset:
        return self._dataset

    def select_variant(self, variant_id: str) -> None:
        """Select a variant and update both panels."""
        self._current_variant = variant_id
        variant_info = self._dataset.get_variant(variant_id)
        if variant_info is None:
            return

        cleavage_data = self._dataset.cleavage_data.get(variant_id, [])

        # Pass modification state if available
        mod_state = None
        if isinstance(self._dataset, EnrichedVariantDataset):
            mod_state = self._dataset.get_modification_state(variant_id)

        self._rna.set_variant(variant_info, cleavage_data, mod_state)
        self._landscape.highlight_variant(variant_id)
        self.variant_selected.emit(variant_id)

    def change_cleavage_site(self, site: int) -> None:
        """Change the displayed cleavage site."""
        self._current_site = site
        self._landscape.update_cleavage_site(self._dataset, site)

        if self._current_variant:
            variant_info = self._dataset.get_variant(self._current_variant)
            if variant_info:
                cleavage_data = self._dataset.cleavage_data.get(
                    self._current_variant, []
                )
                mod_state = None
                if isinstance(self._dataset, EnrichedVariantDataset):
                    mod_state = self._dataset.get_modification_state(self._current_variant)
                self._rna.set_variant(variant_info, cleavage_data, mod_state)

        self.cleavage_site_changed.emit(site)

    def set_dataset(self, dataset: VariantDataset) -> None:
        """Switch to a different dataset (e.g., changing enzyme)."""
        self._dataset = dataset
        self._landscape.build_scatter(dataset, self._current_site)

        if dataset.variants:
            self.select_variant(dataset.variants[0].variant)

        self.enzyme_changed.emit(dataset.enzyme)

    def change_color_mode(self, mode: str) -> None:
        """Change the RNA base coloring mode."""
        self._rna.set_color_mode(mode)
        self.color_mode_changed.emit(mode)

    def change_landscape_mode(self, mode: str) -> None:
        """Switch landscape layout between grid and PCA."""
        self._landscape.set_layout_mode(mode)
        self._landscape.build_scatter(self._dataset, self._current_site)
        if self._current_variant:
            self._landscape.highlight_variant(self._current_variant)
        self.landscape_mode_changed.emit(mode)

    def on_modification_applied(
        self, variant_id: str, position: int, mod_code: str
    ) -> None:
        """Handle a modification being applied — incremental RNA update."""
        if variant_id == self._current_variant and isinstance(
            self._dataset, EnrichedVariantDataset
        ):
            mod_state = self._dataset.get_modification_state(variant_id)
            self._rna.update_modifications_only(mod_state)
        self.modification_applied.emit(variant_id, position, mod_code)
        self.synthesis_updated.emit(variant_id)

    def on_replacement_applied(
        self, variant_id: str, mods: dict[int, str]
    ) -> None:
        """Handle replacement (single/double) — refresh RNA + emit signals."""
        if variant_id == self._current_variant and isinstance(
            self._dataset, EnrichedVariantDataset
        ):
            mod_state = self._dataset.get_modification_state(variant_id)
            self._rna.update_modifications_only(mod_state)
        self.replacement_applied.emit(variant_id, mods)
        self.synthesis_updated.emit(variant_id)
