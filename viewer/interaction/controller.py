"""Cross-panel interaction controller linking RNA view and data landscape."""

from __future__ import annotations

from PySide6.QtCore import QObject, Signal

from viewer.data.schema import VariantDataset
from viewer.landscape.scene import DataLandscapeScene
from viewer.rna3d.scene import RNAHairpinScene


class InteractionController(QObject):
    """Coordinates selection state between the two 3D panels."""

    variant_selected = Signal(str)
    cleavage_site_changed = Signal(int)
    enzyme_changed = Signal(str)

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
        self._rna.set_variant(variant_info, cleavage_data)
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
                self._rna.set_variant(variant_info, cleavage_data)

        self.cleavage_site_changed.emit(site)

    def set_dataset(self, dataset: VariantDataset) -> None:
        """Switch to a different dataset (e.g., changing enzyme)."""
        self._dataset = dataset
        self._landscape.build_scatter(dataset, self._current_site)

        if dataset.variants:
            self.select_variant(dataset.variants[0].variant)

        self.enzyme_changed.emit(dataset.enzyme)
