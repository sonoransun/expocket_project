"""Position x modification type predicted cleavage shift grid."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from matplotlib.figure import Figure

from viewer.encoding.modification_db import MODIFICATIONS_DB

if TYPE_CHECKING:
    from viewer.chemistry.cleavage_predictor import CleavageSitePredictor
    from viewer.chemistry.modification_engine import ModificationEngine

# Cache: variant_id -> impact_matrix
_impact_cache: dict[str, np.ndarray] = {}


def compute_real_impact_matrix(
    variant_id: str,
    predictor: CleavageSitePredictor,
    engine: ModificationEngine,
    seq_length: int,
    mod_codes: list[str],
) -> np.ndarray:
    """Compute actual predicted impact matrix for a variant."""
    if variant_id in _impact_cache:
        return _impact_cache[variant_id]

    matrix = np.zeros((seq_length, len(mod_codes)), dtype=np.float64)
    for j, mod_code in enumerate(mod_codes):
        for pos in range(seq_length):
            applicable = engine.applicable_at_position(variant_id, pos)
            if mod_code not in applicable:
                continue
            ratio = predictor.predict_dc_ratio_shift(variant_id, {pos: mod_code})
            matrix[pos, j] = ratio

    _impact_cache[variant_id] = matrix
    return matrix


def clear_impact_cache() -> None:
    """Clear the impact matrix cache."""
    _impact_cache.clear()


def plot_modification_impact(
    fig: Figure,
    impact_matrix: np.ndarray | None = None,
    seq_length: int = 63,
    mod_codes: list[str] | None = None,
    variant_id: str | None = None,
    predictor: CleavageSitePredictor | None = None,
    engine: ModificationEngine | None = None,
) -> None:
    """Draw a position x modification grid showing predicted cleavage shifts.

    When predictor, engine, and variant_id are provided, computes real
    predicted shifts instead of using demo data.
    """
    fig.clear()

    if mod_codes is None:
        mod_codes = list(MODIFICATIONS_DB.keys())

    n_mods = len(mod_codes)

    # Use real data if engines are available
    if impact_matrix is None and predictor is not None and engine is not None and variant_id is not None:
        impact_matrix = compute_real_impact_matrix(
            variant_id, predictor, engine, seq_length, mod_codes
        )

    if impact_matrix is None:
        # Fallback demo data
        rng = np.random.default_rng(123)
        impact_matrix = rng.normal(0, 0.02, size=(seq_length, n_mods))
        for pos in range(18, min(25, seq_length)):
            impact_matrix[pos, :] *= 3.0
        for j, code in enumerate(mod_codes):
            if code in ("LNA", "2OMe"):
                impact_matrix[:, j] *= 2.0

    impact_matrix = impact_matrix[:seq_length, :n_mods]

    ax = fig.add_axes([0.12, 0.15, 0.78, 0.72])

    vmax = max(abs(impact_matrix.min()), abs(impact_matrix.max()), 0.01)
    im = ax.imshow(impact_matrix.T, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest",
                   origin="lower")

    ax.set_yticks(range(n_mods))
    ax.set_yticklabels(mod_codes, fontsize=6)
    xticks = list(range(0, seq_length, 5))
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=6)
    ax.set_xlabel("RNA position", fontsize=7)
    ax.set_ylabel("Modification", fontsize=7)

    ax.axvspan(19, 23, alpha=0.15, color="red", label="Cleavage zone")
    ax.legend(loc="upper right", fontsize=5, framealpha=0.8)

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Predicted DC ratio shift", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    title = "Modification Impact on Cleavage"
    if variant_id and predictor is not None:
        title += f" — {variant_id}"
    fig.suptitle(title, fontsize=9, y=0.95)
