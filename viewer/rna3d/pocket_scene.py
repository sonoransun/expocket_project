"""Schematic DICER pocket overlay for the RNA 3D view."""

from __future__ import annotations

import numpy as np
import pygfx

from viewer.config import DOMAIN_COLORS
from viewer.encoding.protein_descriptors import DicerPocketModel
from viewer.rna3d.layout import HairpinLayout


class DicerPocketOverlay:
    """Renders a simplified DICER pocket as spheres and contact lines."""

    def __init__(self) -> None:
        self.group = pygfx.Group()
        self._visible = False

    @property
    def visible(self) -> bool:
        return self._visible

    def toggle(self) -> bool:
        self._visible = not self._visible
        self.group.visible = self._visible
        return self._visible

    def set_visible(self, visible: bool) -> None:
        self._visible = visible
        self.group.visible = visible

    def build(self, pocket: DicerPocketModel, layout: HairpinLayout) -> None:
        """Build pocket geometry relative to an RNA layout."""
        self.group.clear()
        if not pocket.residues:
            return

        sphere_geom = pygfx.sphere_geometry(radius=2.0, width_segments=8, height_segments=6)

        for res in pocket.residues:
            # Position residue near the average of its RNA contact positions,
            # offset outward from the backbone
            contact_positions = []
            for rna_pos in res.rna_contact_positions:
                if 0 <= rna_pos < layout.sequence_length:
                    contact_positions.append(layout.base_positions[rna_pos])

            if not contact_positions:
                continue

            center = np.mean(contact_positions, axis=0)
            # Offset outward (away from helix axis) for visual separation
            helix_center = layout.base_positions.mean(axis=0)
            direction = center - helix_center
            norm = np.linalg.norm(direction[:2])
            if norm > 0.1:
                direction[:2] = direction[:2] / norm * 8.0
            else:
                direction[:2] = np.array([8.0, 0.0])
            res_pos = center + direction * 0.5

            # Residue sphere colored by domain
            domain_color = DOMAIN_COLORS.get(res.domain, (0.5, 0.5, 0.5, 0.8))
            material = pygfx.MeshPhongMaterial(
                color=domain_color[:3],
                opacity=domain_color[3] * res.interaction_strength,
            )
            material.transparent = True
            mesh = pygfx.Mesh(sphere_geom, material)
            mesh.local.position = tuple(res_pos)
            scale = 0.5 + res.interaction_strength * 0.8
            mesh.local.scale = (scale, scale, scale)
            self.group.add(mesh)

            # Contact lines from residue to RNA positions
            for rna_pos in res.rna_contact_positions:
                if 0 <= rna_pos < layout.sequence_length:
                    rna_xyz = layout.base_positions[rna_pos].astype(np.float32)
                    res_xyz = np.array(res_pos, dtype=np.float32)
                    positions = np.array([res_xyz, rna_xyz], dtype=np.float32)
                    colors = np.array([
                        [*domain_color[:3], 0.3],
                        [*domain_color[:3], 0.1],
                    ], dtype=np.float32)
                    geom = pygfx.Geometry(positions=positions, colors=colors)
                    line = pygfx.Line(
                        geom, pygfx.LineMaterial(thickness=1.0, color_mode="vertex")
                    )
                    self.group.add(line)

        self.group.visible = self._visible
