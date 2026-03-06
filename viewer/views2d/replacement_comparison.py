"""Grouped bar chart comparing original vs single vs double replacement shifts."""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure


def plot_replacement_comparison(
    fig: Figure,
    original: dict[int, float],
    single_1: dict[int, float] | None = None,
    single_2: dict[int, float] | None = None,
    double: dict[int, float] | None = None,
    labels: list[str] | None = None,
) -> None:
    """Draw grouped bar chart of cleavage accuracy shifts.

    Parameters
    ----------
    original : dict
        Baseline shifts (typically all zeros).
    single_1, single_2, double : dict or None
        Per-site accuracy deltas for single/double replacements.
    labels : list of str or None
        Bar labels. Defaults to ["Original", "Single 1", "Single 2", "Double"].
    """
    fig.clear()
    ax = fig.add_axes([0.12, 0.15, 0.82, 0.72])

    sites = [20, 21, 22, 23]
    n_sites = len(sites)

    groups = [original]
    group_labels = ["Original"]
    colors = ["#888"]

    if single_1 is not None:
        groups.append(single_1)
        group_labels.append("Single 1")
        colors.append("#5AC9A1")

    if single_2 is not None:
        groups.append(single_2)
        group_labels.append("Single 2")
        colors.append("#89CFF0")

    if double is not None:
        groups.append(double)
        group_labels.append("Double")
        colors.append("#ca5cdd")

    if labels is not None:
        group_labels = labels[:len(groups)]

    n_groups = len(groups)
    bar_width = 0.8 / n_groups
    x = np.arange(n_sites)

    for i, (data, label, color) in enumerate(zip(groups, group_labels, colors)):
        values = [data.get(s, 0.0) for s in sites]
        offset = (i - n_groups / 2 + 0.5) * bar_width
        bars = ax.bar(x + offset, values, bar_width, label=label, color=color, alpha=0.85)

        # Delta annotations
        for bar, val in zip(bars, values):
            if abs(val) > 0.001:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.002 * np.sign(val),
                    f"{val:+.3f}", ha="center", va="bottom" if val >= 0 else "top",
                    fontsize=5, color="#ccc",
                )

    ax.set_xticks(x)
    ax.set_xticklabels([f"DC{s}" for s in sites], fontsize=7)
    ax.set_ylabel("Accuracy shift (Δ)", fontsize=7)
    ax.axhline(0, color="#555", linewidth=0.5)
    ax.legend(fontsize=6, framealpha=0.8)

    fig.suptitle("Replacement Comparison", fontsize=9, y=0.95)
