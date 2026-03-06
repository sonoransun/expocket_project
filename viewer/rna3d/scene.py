"""pygfx scene for the RNA hairpin 3D view."""

from __future__ import annotations

import numpy as np
import pygfx

from viewer.config import (
    CLEAVAGE_COLORS_FLOAT,
    CLEAVAGE_SITES,
    MODIFICATION_COLORS,
    NT_COLORS_FLOAT,
    PROPERTY_COLORMAPS,
)
from viewer.data.schema import CleavageRecord, ModificationState, VariantInfo
from viewer.encoding.modification_db import MODIFICATIONS_DB
from viewer.encoding.nucleotide_properties import PROPERTY_NAMES, get_property
from viewer.rna3d.layout import HairpinLayout, compute_hairpin_3d


def _property_to_color(value: float, vmin: float, vmax: float) -> tuple[float, ...]:
    """Map a scalar property value to an RGB color using a blue-red gradient."""
    if vmax == vmin:
        t = 0.5
    else:
        t = (value - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    # Blue -> White -> Red gradient
    if t < 0.5:
        s = t * 2
        return (s, s, 1.0, 1.0)
    else:
        s = (1.0 - t) * 2
        return (1.0, s, s, 1.0)


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
        self._modifications_group = pygfx.Group()
        self.scene.add(self._bases_group)
        self.scene.add(self._backbone_group)
        self.scene.add(self._pairs_group)
        self.scene.add(self._cleavage_group)
        self.scene.add(self._labels_group)
        self.scene.add(self._modifications_group)

        self._layout: HairpinLayout | None = None
        self._color_mode: str = "nucleotide"
        self._current_variant: VariantInfo | None = None
        self._current_mod_state: ModificationState | None = None

    def clear(self) -> None:
        """Remove all visual elements."""
        for group in [
            self._bases_group,
            self._backbone_group,
            self._pairs_group,
            self._cleavage_group,
            self._labels_group,
            self._modifications_group,
        ]:
            group.clear()

    def build_from_variant(
        self,
        variant: VariantInfo,
        cleavage_data: list[CleavageRecord] | None = None,
        mod_state: ModificationState | None = None,
    ) -> None:
        """Build the full 3D scene for a variant."""
        self.clear()
        self._current_variant = variant
        self._current_mod_state = mod_state

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
        if mod_state and mod_state.modifications:
            self._add_modifications(layout, mod_state)

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

        # Precompute property-based colors if needed
        prop_colors = None
        if self._color_mode != "nucleotide" and self._color_mode in PROPERTY_NAMES:
            prop_idx = PROPERTY_NAMES.index(self._color_mode)
            values = [get_property(nt).as_vector()[prop_idx] for nt in layout.base_identities]
            vmin, vmax = min(values), max(values)
            prop_colors = [_property_to_color(v, vmin, vmax) for v in values]

        for i, (pos, nt, region) in enumerate(
            zip(
                layout.base_positions,
                layout.base_identities,
                layout.region_labels,
            )
        ):
            if prop_colors is not None:
                color = prop_colors[i]
            else:
                color = NT_COLORS_FLOAT.get(nt.upper(), (0.5, 0.5, 0.5, 1.0))

            # Highlight randomized nucleotides (last 3 bases) with larger radius
            is_randomized = i >= layout.sequence_length - n_randomized
            radius_scale = 1.4 if is_randomized else 1.0

            # Highlight loop region slightly brighter
            if region == "loop" and prop_colors is None:
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

    def _add_modifications(
        self, layout: HairpinLayout, mod_state: ModificationState
    ) -> None:
        """Render modification overlays at modified positions."""
        for pos_idx, mod_code in mod_state.modifications.items():
            if pos_idx < 0 or pos_idx >= layout.sequence_length:
                continue
            mod = MODIFICATIONS_DB.get(mod_code)
            if mod is None:
                continue

            base_pos = layout.base_positions[pos_idx]
            color = MODIFICATION_COLORS.get(mod_code, (0.8, 0.8, 0.8, 0.85))

            # Different geometry per modification type
            shape = mod.display_shape
            if shape == "cube":
                geom = pygfx.box_geometry(3.5, 3.5, 3.5)
            elif shape == "octahedron":
                geom = pygfx.sphere_geometry(radius=2.2, width_segments=4, height_segments=3)
            elif shape == "dodecahedron":
                geom = pygfx.sphere_geometry(radius=2.5, width_segments=6, height_segments=5)
            else:  # cone / default
                geom = pygfx.cone_geometry(radius=2.0, height=3.5, radial_segments=12)

            material = pygfx.MeshPhongMaterial(
                color=color[:3],
                opacity=color[3],
                emissive=tuple(c * 0.3 for c in color[:3]),
            )
            material.transparent = True
            mesh = pygfx.Mesh(geom, material)
            mesh.local.position = tuple(base_pos)
            self._modifications_group.add(mesh)

    def set_color_mode(self, mode: str) -> None:
        """Change base coloring mode and rebuild bases only."""
        self._color_mode = mode
        if self._layout is not None and self._current_variant is not None:
            self._bases_group.clear()
            self._add_bases(self._layout, self._current_variant)

    def set_variant(
        self,
        variant: VariantInfo,
        cleavage_data: list[CleavageRecord] | None = None,
        mod_state: ModificationState | None = None,
    ) -> None:
        """Update the scene for a newly selected variant."""
        self.build_from_variant(variant, cleavage_data, mod_state)
