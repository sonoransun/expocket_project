"""Position x modification type predicted cleavage shift grid."""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from viewer.encoding.modification_db import MODIFICATIONS_DB


def plot_modification_impact(
    fig: Figure,
    impact_matrix: np.ndarray | None = None,
    seq_length: int = 63,
    mod_codes: list[str] | None = None,
) -> None:
    """Draw a position x modification grid showing predicted cleavage shifts.

    Parameters
    ----------
    impact_matrix : (n_positions, n_modifications) array or None
        Predicted DC21/DC22 ratio shift. If None, generates a demo matrix.
    seq_length : int
        Length of the RNA sequence.
    mod_codes : list of str or None
        Modification codes (column labels). Defaults to all in DB.
    """
    fig.clear()

    if mod_codes is None:
        mod_codes = list(MODIFICATIONS_DB.keys())

    n_mods = len(mod_codes)

    if impact_matrix is None:
        # Generate a plausible demo: modifications near cleavage sites have larger effects
        rng = np.random.default_rng(123)
        impact_matrix = rng.normal(0, 0.02, size=(seq_length, n_mods))
        # Amplify near cleavage zone (positions 18-24)
        for pos in range(18, min(25, seq_length)):
            impact_matrix[pos, :] *= 3.0
        # LNA and 2OMe have strongest effects
        for j, code in enumerate(mod_codes):
            if code in ("LNA", "2OMe"):
                impact_matrix[:, j] *= 2.0

    # Trim to actual size
    impact_matrix = impact_matrix[:seq_length, :n_mods]

    ax = fig.add_axes([0.12, 0.15, 0.78, 0.72])

    vmax = max(abs(impact_matrix.min()), abs(impact_matrix.max()), 0.1)
    im = ax.imshow(impact_matrix.T, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest",
                   origin="lower")

    # Labels
    ax.set_yticks(range(n_mods))
    ax.set_yticklabels(mod_codes, fontsize=6)
    xticks = list(range(0, seq_length, 5))
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=6)
    ax.set_xlabel("RNA position", fontsize=7)
    ax.set_ylabel("Modification", fontsize=7)

    # Cleavage zone highlight
    ax.axvspan(19, 23, alpha=0.15, color="red", label="Cleavage zone")
    ax.legend(loc="upper right", fontsize=5, framealpha=0.8)

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Predicted DC ratio shift", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle("Modification Impact on Cleavage", fontsize=9, y=0.95)
