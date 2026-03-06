"""pygfx scene for the RNA hairpin 3D view."""

from __future__ import annotations

import numpy as np
import pygfx

from viewer.config import (
    CLEAVAGE_COLORS_FLOAT,
    CLEAVAGE_SITES,
    NT_COLORS_FLOAT,
)
from viewer.data.schema import CleavageRecord, VariantInfo
from viewer.rna3d.layout import HairpinLayout, compute_hairpin_3d


class RNAHairpinScene:
    """Manages the pygfx scene graph for the RNA structure panel."""

    def __init__(self) -> None:
        self.scene = pygfx.Scene()
        self.camera = pygfx.PerspectiveCamera(fov=50)
        self.camera.local.position = (0, -50, 100)

        # Lighting
        self.scene.add(pygfx.AmbientLight(intensity=0.4))
        light = pygfx.DirectionalLight(intensity=0.8)
        light.local.position = (50, 80, 100)
        self.scene.add(light)

        # Background
        self.scene.add(
            pygfx.Background.from_color("#1a1a2e", "#16213e")
        )

        # Groups for different element types
        self._bases_group = pygfx.Group()
        self._backbone_group = pygfx.Group()
        self._pairs_group = pygfx.Group()
        self._cleavage_group = pygfx.Group()
        self._labels_group = pygfx.Group()
        self.scene.add(self._bases_group)
        self.scene.add(self._backbone_group)
        self.scene.add(self._pairs_group)
        self.scene.add(self._cleavage_group)
        self.scene.add(self._labels_group)

        self._layout: HairpinLayout | None = None

    def clear(self) -> None:
        """Remove all visual elements."""
        for group in [
            self._bases_group,
            self._backbone_group,
            self._pairs_group,
            self._cleavage_group,
            self._labels_group,
        ]:
            group.clear()

    def build_from_variant(
        self,
        variant: VariantInfo,
        cleavage_data: list[CleavageRecord] | None = None,
    ) -> None:
        """Build the full 3D scene for a variant."""
        self.clear()

        layout = compute_hairpin_3d(
            variant.pre_mirna_sequence,
            variant.concrete_struct,
            variant.flanking_length_5p,
        )
        self._layout = layout

        self._add_backbone(layout)
        self._add_bases(layout, variant)
        self._add_base_pairs(layout)
        self._add_cleavage_markers(layout, cleavage_data)

        # Center camera on the structure
        center = layout.base_positions.mean(axis=0)
        extent = np.ptp(layout.base_positions, axis=0).max()
        self.camera.local.position = (
            center[0],
            center[1] - extent * 0.3,
            center[2] + extent * 1.5,
        )
        self.camera.show_object(self.scene, view_dir=(0, 0.2, -1))

    def _add_backbone(self, layout: HairpinLayout) -> None:
        """Render the backbone as a line strip."""
        positions = layout.backbone_spline.astype(np.float32)
        colors = np.full((len(positions), 4), (0.6, 0.6, 0.7, 0.8), dtype=np.float32)

        geom = pygfx.Geometry(positions=positions, colors=colors)
        material = pygfx.LineMaterial(thickness=3.0, color_mode="vertex")
        line = pygfx.Line(geom, material)
        self._backbone_group.add(line)

    def _add_bases(self, layout: HairpinLayout, variant: VariantInfo) -> None:
        """Render each nucleotide base as a colored sphere."""
        sphere_geom = pygfx.sphere_geometry(radius=1.8, width_segments=12, height_segments=8)
        n_randomized = len(variant.randomized_nts)

        for i, (pos, nt, region) in enumerate(
            zip(
                layout.base_positions,
                layout.base_identities,
                layout.region_labels,
            )
        ):
            color = NT_COLORS_FLOAT.get(nt.upper(), (0.5, 0.5, 0.5, 1.0))

            # Highlight randomized nucleotides (last 3 bases) with larger radius
            is_randomized = i >= layout.sequence_length - n_randomized
            radius_scale = 1.4 if is_randomized else 1.0

            # Highlight loop region slightly brighter
            if region == "loop":
                color = tuple(min(1.0, c * 1.2) for c in color[:3]) + (color[3],)

            material = pygfx.MeshPhongMaterial(color=color[:3])
            mesh = pygfx.Mesh(sphere_geom, material)
            mesh.local.position = tuple(pos)
            mesh.local.scale = (radius_scale, radius_scale, radius_scale)
            self._bases_group.add(mesh)

    def _add_base_pairs(self, layout: HairpinLayout) -> None:
        """Draw thin lines between paired bases (hydrogen bonds)."""
        for i5, i3 in layout.pair_indices:
            pos5 = layout.base_positions[i5].astype(np.float32)
            pos3 = layout.base_positions[i3].astype(np.float32)
            positions = np.array([pos5, pos3], dtype=np.float32)
            colors = np.full((2, 4), (0.4, 0.4, 0.5, 0.5), dtype=np.float32)
            geom = pygfx.Geometry(positions=positions, colors=colors)
            material = pygfx.LineMaterial(thickness=1.0, color_mode="vertex")
            line = pygfx.Line(geom, material)
            self._pairs_group.add(line)

    def _add_cleavage_markers(
        self,
        layout: HairpinLayout,
        cleavage_data: list[CleavageRecord] | None = None,
    ) -> None:
        """Add semi-transparent colored discs at cleavage sites."""
        accuracy_map: dict[int, float] = {}
        if cleavage_data:
            for rec in cleavage_data:
                accuracy_map[rec.cleavage_site] = rec.mean_accuracy

        for site in CLEAVAGE_SITES:
            if site not in layout.cleavage_site_positions:
                continue
            pos = layout.cleavage_site_positions[site]
            color = CLEAVAGE_COLORS_FLOAT.get(site, (0.5, 0.5, 0.5, 0.4))

            # Scale disc by accuracy if available
            acc = accuracy_map.get(site, 0.3)
            disc_radius = 4.0 + acc * 8.0

            geom = pygfx.cylinder_geometry(
                radius_top=disc_radius,
                radius_bottom=disc_radius,
                height=0.5,
                radial_segments=32,
            )
            material = pygfx.MeshPhongMaterial(
                color=color[:3],
                opacity=color[3] * (0.3 + acc * 0.7),
            )
            material.transparent = True
            mesh = pygfx.Mesh(geom, material)
            mesh.local.position = tuple(pos)
            self._cleavage_group.add(mesh)

    def set_variant(
        self,
        variant: VariantInfo,
        cleavage_data: list[CleavageRecord] | None = None,
    ) -> None:
        """Update the scene for a newly selected variant."""
        self.build_from_variant(variant, cleavage_data)
