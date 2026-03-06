"""pygfx scene for the 3D data landscape of 256 variants."""

from __future__ import annotations

import numpy as np
import pygfx

from viewer.config import CLEAVAGE_COLORS_FLOAT, GROUPS, NT_COLORS_FLOAT
from viewer.data.schema import EnrichedVariantDataset, VariantDataset


# Grid layout: 4 groups on X, 64 variants per group spread on Y
_GROUP_X = {g: i * 20.0 for i, g in enumerate(GROUPS)}


class DataLandscapeScene:
    """Manages the pygfx scene graph for the variant data landscape."""

    def __init__(self) -> None:
        self.scene = pygfx.Scene()
        self.camera = pygfx.PerspectiveCamera(fov=50)

        # Lighting
        self.scene.add(pygfx.AmbientLight(intensity=0.5))
        light = pygfx.DirectionalLight(intensity=0.7)
        light.local.position = (60, 80, 100)
        self.scene.add(light)

        self.scene.add(
            pygfx.Background.from_color("#1a1a2e", "#0f3460")
        )

        self._points_group = pygfx.Group()
        self._axes_group = pygfx.Group()
        self._highlight_group = pygfx.Group()
        self.scene.add(self._points_group)
        self.scene.add(self._axes_group)
        self.scene.add(self._highlight_group)

        self._variant_meshes: dict[str, pygfx.Mesh] = {}
        self._variant_positions: dict[str, np.ndarray] = {}
        self._mesh_to_variant: dict[int, str] = {}  # id(mesh) -> variant_id
        self._dataset: VariantDataset | None = None
        self._layout_mode: str = "grid"  # "grid" or "pca"

    def clear(self) -> None:
        self._points_group.clear()
        self._axes_group.clear()
        self._highlight_group.clear()
        self._variant_meshes.clear()
        self._variant_positions.clear()
        self._mesh_to_variant.clear()

    def set_layout_mode(self, mode: str) -> None:
        """Switch between 'grid' and 'pca' layout modes."""
        if mode in ("grid", "pca"):
            self._layout_mode = mode

    def build_scatter(
        self, dataset: VariantDataset, cleavage_site: int = 21
    ) -> None:
        """Build a 3D scatter plot for all variants at a given cleavage site."""
        self.clear()
        self._dataset = dataset

        if self._layout_mode == "pca" and isinstance(dataset, EnrichedVariantDataset):
            self._build_pca_scatter(dataset, cleavage_site)
        else:
            self._build_grid_scatter(dataset, cleavage_site)

        self._add_axes(cleavage_site)
        self._center_camera()

    def _build_pca_scatter(
        self, dataset: EnrichedVariantDataset, cleavage_site: int
    ) -> None:
        """Build scatter using PCA of summary features for X/Y, accuracy for Z."""
        sphere_geom = pygfx.sphere_geometry(
            radius=1.2, width_segments=10, height_segments=6
        )

        features = dataset.summary_features
        if features is None or features.shape[0] == 0:
            return

        # PCA via SVD (no sklearn needed)
        centered = features - features.mean(axis=0)
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        pc = centered @ Vt[:2].T  # (N, 2) — first two PCs

        # Scale PCs to reasonable visual range
        pc_range = np.ptp(pc, axis=0)
        pc_range[pc_range == 0] = 1.0
        pc_scaled = pc / pc_range * 40.0

        for i, variant in enumerate(dataset.variants):
            vid = variant.variant
            group = variant.group

            x = float(pc_scaled[i, 0])
            y = float(pc_scaled[i, 1])

            rec = dataset.get_cleavage(vid, cleavage_site)
            z = (rec.mean_accuracy if rec else 0.0) * 40.0

            color = NT_COLORS_FLOAT.get(group, (0.5, 0.5, 0.5, 1.0))
            acc = rec.mean_accuracy if rec else 0.0
            bright = 0.5 + acc * 0.5
            color_mod = (color[0] * bright, color[1] * bright, color[2] * bright)

            material = pygfx.MeshPhongMaterial(color=color_mod)
            mesh = pygfx.Mesh(sphere_geom, material)
            pos = np.array([x, y, z])
            mesh.local.position = tuple(pos)
            self._points_group.add(mesh)
            self._variant_meshes[vid] = mesh
            self._variant_positions[vid] = pos
            self._mesh_to_variant[id(mesh)] = vid

    def _build_grid_scatter(
        self, dataset: VariantDataset, cleavage_site: int
    ) -> None:
        """Build scatter using the original grid layout."""
        sphere_geom = pygfx.sphere_geometry(
            radius=1.2, width_segments=10, height_segments=6
        )

        group_counts: dict[str, int] = {g: 0 for g in GROUPS}

        for variant in dataset.variants:
            vid = variant.variant
            group = variant.group
            idx_in_group = group_counts.get(group, 0)
            group_counts[group] = idx_in_group + 1

            row = idx_in_group // 8
            col = idx_in_group % 8

            x = _GROUP_X.get(group, 0) + col * 2.0
            y = row * 2.0

            rec = dataset.get_cleavage(vid, cleavage_site)
            z = (rec.mean_accuracy if rec else 0.0) * 40.0

            color = NT_COLORS_FLOAT.get(group, (0.5, 0.5, 0.5, 1.0))
            acc = rec.mean_accuracy if rec else 0.0
            bright = 0.5 + acc * 0.5
            color_mod = (
                color[0] * bright,
                color[1] * bright,
                color[2] * bright,
            )

            material = pygfx.MeshPhongMaterial(color=color_mod)
            mesh = pygfx.Mesh(sphere_geom, material)
            pos = np.array([x, y, z])
            mesh.local.position = tuple(pos)
            self._points_group.add(mesh)
            self._variant_meshes[vid] = mesh
            self._variant_positions[vid] = pos
            self._mesh_to_variant[id(mesh)] = vid

    def _add_axes(self, cleavage_site: int) -> None:
        """Add axis lines and group labels."""
        # X axis line
        x_line = np.array([[0, 0, 0], [80, 0, 0]], dtype=np.float32)
        geom = pygfx.Geometry(
            positions=x_line,
            colors=np.full((2, 4), (0.5, 0.5, 0.5, 0.6), dtype=np.float32),
        )
        self._axes_group.add(
            pygfx.Line(geom, pygfx.LineMaterial(thickness=1.5, color_mode="vertex"))
        )

        # Y axis line
        y_line = np.array([[0, 0, 0], [0, 16, 0]], dtype=np.float32)
        geom = pygfx.Geometry(
            positions=y_line,
            colors=np.full((2, 4), (0.5, 0.5, 0.5, 0.6), dtype=np.float32),
        )
        self._axes_group.add(
            pygfx.Line(geom, pygfx.LineMaterial(thickness=1.5, color_mode="vertex"))
        )

        # Z axis line
        z_line = np.array([[0, 0, 0], [0, 0, 45]], dtype=np.float32)
        geom = pygfx.Geometry(
            positions=z_line,
            colors=np.full((2, 4), (0.5, 0.5, 0.5, 0.6), dtype=np.float32),
        )
        self._axes_group.add(
            pygfx.Line(geom, pygfx.LineMaterial(thickness=1.5, color_mode="vertex"))
        )

        # Ground grid
        grid_lines = []
        for g_idx, group in enumerate(GROUPS):
            gx = _GROUP_X[group]
            for col in range(8):
                x = gx + col * 2.0
                grid_lines.append([x, 0, 0])
                grid_lines.append([x, 14, 0])
        if grid_lines:
            grid_pos = np.array(grid_lines, dtype=np.float32)
            grid_colors = np.full(
                (len(grid_lines), 4), (0.2, 0.2, 0.3, 0.3), dtype=np.float32
            )
            geom = pygfx.Geometry(positions=grid_pos, colors=grid_colors)
            self._axes_group.add(
                pygfx.Line(
                    geom,
                    pygfx.LineSegmentMaterial(thickness=0.5, color_mode="vertex"),
                )
            )

    def _center_camera(self) -> None:
        if not self._variant_positions:
            return
        all_pos = np.array(list(self._variant_positions.values()))
        center = all_pos.mean(axis=0)
        extent = np.ptp(all_pos, axis=0).max()
        self.camera.local.position = (
            center[0],
            center[1] - extent * 0.5,
            center[2] + extent * 1.2,
        )
        self.camera.show_object(self.scene, view_dir=(0, 0.3, -1))

    def highlight_variant(self, variant_id: str) -> None:
        """Highlight a specific variant with a ring/glow."""
        self._highlight_group.clear()
        pos = self._variant_positions.get(variant_id)
        if pos is None:
            return

        highlight_geom = pygfx.sphere_geometry(
            radius=2.5, width_segments=16, height_segments=12
        )
        material = pygfx.MeshPhongMaterial(
            color=(1.0, 1.0, 0.3),
            emissive=(0.8, 0.8, 0.0),
            opacity=0.4,
        )
        material.transparent = True
        highlight = pygfx.Mesh(highlight_geom, material)
        highlight.local.position = tuple(pos)
        self._highlight_group.add(highlight)

    def update_cleavage_site(
        self, dataset: VariantDataset, site: int
    ) -> None:
        """Rebuild the scatter for a different cleavage site."""
        self.build_scatter(dataset, site)

    def get_variant_for_mesh(self, mesh: pygfx.Mesh) -> str | None:
        """Look up variant_id by mesh object (for pick events)."""
        return self._mesh_to_variant.get(id(mesh))

    def get_nearest_variant(self, world_pos: np.ndarray) -> str | None:
        """Find the variant closest to a world position (for picking)."""
        if not self._variant_positions:
            return None
        best_dist = float("inf")
        best_id = None
        for vid, pos in self._variant_positions.items():
            d = float(np.linalg.norm(pos - world_pos))
            if d < best_dist:
                best_dist = d
                best_id = vid
        return best_id if best_dist < 5.0 else None
