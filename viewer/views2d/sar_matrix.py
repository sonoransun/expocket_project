"""Feature-vs-cleavage correlation matrix for SAR analysis."""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from viewer.data.schema import EnrichedVariantDataset
from viewer.encoding.property_calculator import SUMMARY_FEATURE_NAMES


def plot_sar_matrix(
    fig: Figure,
    dataset: EnrichedVariantDataset,
    cleavage_site: int = 21,
) -> None:
    """Draw a correlation matrix: features vs cleavage metrics."""
    fig.clear()

    features = dataset.summary_features
    if features is None or features.shape[0] == 0:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return

    # Gather cleavage accuracy vector
    n = len(dataset.variants)
    accuracy = np.zeros(n)
    for i, v in enumerate(dataset.variants):
        rec = dataset.get_cleavage(v.variant, cleavage_site)
        if rec is not None:
            accuracy[i] = rec.mean_accuracy

    # Compute Pearson correlation between each feature and accuracy
    n_feat = features.shape[1]
    correlations = np.zeros(n_feat)
    p_values = np.zeros(n_feat)

    for j in range(n_feat):
        x = features[:, j]
        if x.std() == 0 or accuracy.std() == 0:
            continue
        r = np.corrcoef(x, accuracy)[0, 1]
        correlations[j] = r
        # t-statistic for significance
        t = r * np.sqrt((n - 2) / (1 - r * r + 1e-12))
        # Approximate p-value using normal for large n
        p_values[j] = 2 * np.exp(-0.5 * t * t) / np.sqrt(2 * np.pi) if abs(t) < 10 else 0.0

    # Reshape into 4 groups x 12 properties for display
    feat_names = SUMMARY_FEATURE_NAMES[:n_feat]
    n_props = 12
    n_groups = n_feat // n_props
    corr_matrix = correlations[:n_groups * n_props].reshape(n_groups, n_props)
    pval_matrix = p_values[:n_groups * n_props].reshape(n_groups, n_props)

    ax = fig.add_axes([0.18, 0.15, 0.72, 0.72])

    im = ax.imshow(corr_matrix, aspect="auto", cmap="RdBu_r",
                   vmin=-0.5, vmax=0.5, interpolation="nearest")

    # Significance annotations
    for i in range(n_groups):
        for j in range(n_props):
            stars = ""
            if pval_matrix[i, j] < 0.001:
                stars = "***"
            elif pval_matrix[i, j] < 0.01:
                stars = "**"
            elif pval_matrix[i, j] < 0.05:
                stars = "*"
            if stars:
                ax.text(j, i, stars, ha="center", va="center", fontsize=6, color="k")

    # Labels
    prop_labels = [n.split("_", 1)[-1] if "_" in n else n for n in feat_names[:n_props]]
    group_labels = ["mean", "5'", "3' rand", "clv zone"][:n_groups]

    ax.set_xticks(range(n_props))
    ax.set_xticklabels(prop_labels, rotation=90, fontsize=6)
    ax.set_yticks(range(n_groups))
    ax.set_yticklabels(group_labels, fontsize=7)

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Pearson r", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle(f"SAR Correlation — DC{cleavage_site} Accuracy", fontsize=9, y=0.95)
