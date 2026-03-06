"""DICER pocket residue x RNA position contact map."""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from viewer.config import CLEAVAGE_SITES, DOMAIN_COLORS
from viewer.encoding.protein_descriptors import DicerPocketModel


def plot_contact_map(
    fig: Figure,
    pocket: DicerPocketModel | None,
    seq_length: int = 63,
) -> None:
    """Draw a contact map: DICER residues (Y) vs RNA positions (X)."""
    fig.clear()

    if pocket is None or not pocket.residues:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No pocket data", ha="center", va="center")
        return

    mat = pocket.contact_matrix(seq_length)
    n_res = len(pocket.residues)

    ax = fig.add_axes([0.18, 0.12, 0.72, 0.78])

    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=1.0, interpolation="nearest",
                   origin="lower")

    # Domain color bar on the left
    domains = pocket.domains
    domain_to_idx = {d: i for i, d in enumerate(domains)}
    for i, res in enumerate(pocket.residues):
        color = DOMAIN_COLORS.get(res.domain, (0.5, 0.5, 0.5, 0.8))
        ax.barh(i, -2, height=0.8, left=-3, color=color[:3], clip_on=False)

    # Residue labels
    res_labels = [f"{r.amino_acid}{r.residue_id}" for r in pocket.residues]
    ax.set_yticks(range(n_res))
    ax.set_yticklabels(res_labels, fontsize=5)

    # RNA position labels (show every 5th)
    xticks = list(range(0, seq_length, 5))
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=6)
    ax.set_xlabel("RNA position", fontsize=7)
    ax.set_ylabel("DICER residue", fontsize=7)

    # Mark cleavage sites
    for site in CLEAVAGE_SITES:
        if site < seq_length:
            ax.axvline(site, color="red", linewidth=0.5, linestyle="--", alpha=0.6)

    # Domain legend
    for domain in domains:
        color = DOMAIN_COLORS.get(domain, (0.5, 0.5, 0.5, 0.8))
        ax.plot([], [], "s", color=color[:3], label=domain, markersize=5)
    ax.legend(loc="upper right", fontsize=5, framealpha=0.8)

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Interaction strength", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle(f"DICER Pocket Contact Map ({pocket.enzyme})", fontsize=9, y=0.97)
