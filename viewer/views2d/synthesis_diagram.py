"""Synthesis pathway visualization: sequence bead diagram and yield curve."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from matplotlib.figure import Figure

if TYPE_CHECKING:
    from viewer.data.schema import SynthesisPlan

from viewer.config import MODIFICATION_COLORS, NT_COLORS


def plot_synthesis_diagram(fig: Figure, plan: SynthesisPlan | None) -> None:
    """Draw synthesis as bead-on-string diagram (top) and yield curve (bottom)."""
    fig.clear()

    if plan is None or not plan.steps:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No synthesis plan — apply modifications first",
                ha="center", va="center", color="#ccc", fontsize=9)
        return

    seq = plan.sequence
    n = len(seq)

    ax_beads = fig.add_axes([0.05, 0.55, 0.90, 0.38])
    ax_yield = fig.add_axes([0.10, 0.10, 0.82, 0.38])

    # --- Bead diagram (top) ---
    # Show nucleotides as colored circles along a line
    x = np.arange(n)
    y = np.zeros(n)

    # Backbone line
    ax_beads.plot(x, y, "-", color="#555", linewidth=1, zorder=1)

    # Draw each nucleotide as a circle
    for i, nt in enumerate(seq):
        color = NT_COLORS.get(nt.upper(), "#888")
        size = 80
        mod = plan.modifications.get(i)

        ax_beads.scatter(i, 0, s=size, c=color, edgecolors="#333",
                         linewidths=0.5, zorder=2)

        # Modification ring
        if mod:
            mod_color_rgba = MODIFICATION_COLORS.get(mod, (0.8, 0.8, 0.8, 0.85))
            mod_color = mod_color_rgba[:3]
            ax_beads.scatter(i, 0, s=200, facecolors="none",
                             edgecolors=mod_color, linewidths=2, zorder=3)
            ax_beads.text(i, 0.6, mod, ha="center", va="bottom",
                          fontsize=4, color=mod_color, rotation=45)

    # Labels
    ax_beads.set_xlim(-1, n)
    ax_beads.set_ylim(-1.2, 1.5)
    ax_beads.set_yticks([])
    # Show position ticks every 10
    ticks = list(range(0, n, 10))
    ax_beads.set_xticks(ticks)
    ax_beads.set_xticklabels([str(t) for t in ticks], fontsize=5)
    ax_beads.text(0, -0.8, "5'", fontsize=7, ha="center", color="#aaa")
    ax_beads.text(n - 1, -0.8, "3'", fontsize=7, ha="center", color="#aaa")
    ax_beads.set_title(
        f"Sequence: {n} nt, {len(plan.modifications)} modifications",
        fontsize=7, color="#ccc",
    )

    # Cleavage zone marker
    ax_beads.axvspan(19, 23, alpha=0.1, color="red")

    # --- Yield curve (bottom) ---
    # Steps are in synthesis order (3'→5'), x-axis = step number
    step_nums = list(range(1, len(plan.steps) + 1))
    yields = [s.cumulative_yield * 100 for s in plan.steps]

    ax_yield.plot(step_nums, yields, "-", color="#5AC9A1", linewidth=1.5)
    ax_yield.fill_between(step_nums, yields, alpha=0.15, color="#5AC9A1")

    # Mark modification steps
    for i, step in enumerate(plan.steps):
        if step.modification:
            ax_yield.plot(i + 1, step.cumulative_yield * 100, "o",
                          color=MODIFICATION_COLORS.get(step.modification, (0.8, 0.8, 0.8, 0.85))[:3],
                          markersize=5, zorder=3)

    # Final yield annotation
    final_yield = plan.total_yield * 100
    ax_yield.annotate(
        f"{final_yield:.1f}%",
        xy=(len(plan.steps), final_yield),
        fontsize=8, color="#5AC9A1", fontweight="bold",
        xytext=(-30, 10), textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color="#5AC9A1", lw=0.8),
    )

    ax_yield.set_xlabel("Synthesis step (3'→5')", fontsize=7)
    ax_yield.set_ylabel("Cumulative yield (%)", fontsize=7)
    ax_yield.set_ylim(0, 105)
    ax_yield.set_xlim(1, len(plan.steps))

    # Cost annotation
    ax_yield.text(
        0.02, 0.95,
        f"Cost: {plan.total_cost_factor:.0f}× | Scale: {plan.scale_nmol} nmol",
        transform=ax_yield.transAxes, fontsize=6, color="#aaa", va="top",
    )

    if plan.incompatibilities:
        warn_text = " | ".join(plan.incompatibilities[:2])
        if len(warn_text) > 60:
            warn_text = warn_text[:57] + "..."
        ax_yield.text(
            0.02, 0.85, f"⚠ {warn_text}",
            transform=ax_yield.transAxes, fontsize=5, color="#f44",
            va="top",
        )

    fig.suptitle(f"Synthesis Plan — {plan.variant_id}", fontsize=9, y=0.98)
