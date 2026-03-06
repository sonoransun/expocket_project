"""Clustered heatmap of physicochemical properties across variants."""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from viewer.data.schema import EnrichedVariantDataset
from viewer.encoding.property_calculator import SUMMARY_FEATURE_NAMES


def plot_property_heatmap(
    fig: Figure,
    dataset: EnrichedVariantDataset,
    max_features: int = 48,
) -> None:
    """Draw a variant x feature clustered heatmap on *fig*.

    Uses Ward linkage on z-score normalized summary features.
    """
    fig.clear()

    features = dataset.summary_features
    if features is None or features.shape[0] == 0:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return

    n_variants, n_feat = features.shape
    n_feat = min(n_feat, max_features)
    data = features[:, :n_feat].copy()

    # Z-score normalize columns
    col_std = data.std(axis=0)
    col_std[col_std == 0] = 1.0
    z_data = (data - data.mean(axis=0)) / col_std

    # Hierarchical clustering of rows (variants)
    row_dist = pdist(z_data, metric="euclidean")
    row_link = linkage(row_dist, method="ward")

    # Layout: dendrogram on left, heatmap on right
    ax_dendro = fig.add_axes([0.02, 0.05, 0.12, 0.88])
    ax_heat = fig.add_axes([0.16, 0.05, 0.78, 0.88])

    dn = dendrogram(row_link, orientation="left", ax=ax_dendro, no_labels=True,
                     color_threshold=0, above_threshold_color="0.5")
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])
    ax_dendro.spines[:].set_visible(False)

    # Reorder rows by dendrogram
    row_order = dn["leaves"]
    ordered = z_data[row_order]

    im = ax_heat.imshow(ordered, aspect="auto", cmap="RdBu_r",
                        vmin=-2.5, vmax=2.5, interpolation="nearest")

    # Feature labels on x-axis
    feat_names = SUMMARY_FEATURE_NAMES[:n_feat]
    ax_heat.set_xticks(range(n_feat))
    ax_heat.set_xticklabels(feat_names, rotation=90, fontsize=5)
    ax_heat.set_yticks([])
    ax_heat.set_ylabel(f"{n_variants} variants (clustered)")

    cb = fig.colorbar(im, ax=ax_heat, fraction=0.02, pad=0.01)
    cb.set_label("z-score", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle("Physicochemical Property Heatmap", fontsize=9, y=0.98)
