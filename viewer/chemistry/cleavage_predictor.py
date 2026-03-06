"""Ridge regression predictor: property perturbations -> cleavage shift."""

from __future__ import annotations

import numpy as np

from viewer.data.schema import EnrichedVariantDataset
from viewer.encoding.property_calculator import (
    compute_modified_properties,
    compute_summary_features,
)


class CleavageSitePredictor:
    """Predicts cleavage accuracy shifts from physicochemical features.

    Uses ridge regression (closed-form) fit on the dataset's summary features.
    """

    def __init__(self, dataset: EnrichedVariantDataset, alpha: float = 1.0) -> None:
        self._dataset = dataset
        self._alpha = alpha
        self._weights: dict[int, np.ndarray] = {}  # site -> (F,) coefficients
        self._intercepts: dict[int, float] = {}
        self._fit()

    def _fit(self) -> None:
        """Fit ridge regression for each cleavage site."""
        features = self._dataset.summary_features
        if features is None or features.shape[0] == 0:
            return

        X = features.copy()
        n, f = X.shape

        # Center and scale
        self._mean = X.mean(axis=0)
        self._std = X.std(axis=0)
        self._std[self._std == 0] = 1.0
        X_norm = (X - self._mean) / self._std

        for site in [20, 21, 22, 23]:
            y = np.zeros(n)
            for i, v in enumerate(self._dataset.variants):
                rec = self._dataset.get_cleavage(v.variant, site)
                if rec is not None:
                    y[i] = rec.mean_accuracy

            y_mean = y.mean()
            y_c = y - y_mean

            # Ridge: w = (X'X + alpha*I)^-1 X'y
            XtX = X_norm.T @ X_norm + self._alpha * np.eye(f)
            Xty = X_norm.T @ y_c
            w = np.linalg.solve(XtX, Xty)

            self._weights[site] = w
            self._intercepts[site] = y_mean

    def predict_shift(
        self,
        variant_id: str,
        modifications: dict[int, str],
    ) -> dict[int, float]:
        """Predict accuracy change at each cleavage site from modifications.

        Returns {site: delta_accuracy} for sites 20-23.
        """
        variant = self._dataset.get_variant(variant_id)
        if variant is None or not self._weights:
            return {s: 0.0 for s in [20, 21, 22, 23]}

        # Baseline features
        from viewer.encoding.property_calculator import compute_summary_features as _csf
        from viewer.data.schema import VariantDataset

        baseline_ds = VariantDataset(variants=[variant])
        baseline_feat = _csf(baseline_ds)[0]

        # Modified features — create a temporary single-variant dataset
        mod_props = compute_modified_properties(
            variant.pre_mirna_sequence, modifications
        )
        # Recompute summary from modified property matrix
        mod_feat = self._summarize_single(mod_props)

        result = {}
        for site in [20, 21, 22, 23]:
            w = self._weights.get(site)
            if w is None:
                result[site] = 0.0
                continue
            base_pred = float(((baseline_feat - self._mean) / self._std) @ w)
            mod_pred = float(((mod_feat - self._mean) / self._std) @ w)
            result[site] = mod_pred - base_pred
        return result

    def predict_dc_ratio_shift(
        self,
        variant_id: str,
        modifications: dict[int, str],
    ) -> float:
        """Predict DC21/DC22 ratio change from modifications."""
        shifts = self.predict_shift(variant_id, modifications)
        return shifts.get(21, 0.0) - shifts.get(22, 0.0)

    @staticmethod
    def _summarize_single(props: np.ndarray) -> np.ndarray:
        """Compute 48-feature summary from a single (N, 12) matrix."""
        result = np.zeros(48, dtype=np.float64)
        n = props.shape[0]
        result[0:12] = props.mean(axis=0)
        if n > 0:
            result[12:24] = props[0]
        if n >= 3:
            result[24:36] = props[-3:].mean(axis=0)
        clv_start = min(19, n - 1)
        clv_end = min(23, n)
        if clv_end > clv_start:
            result[36:48] = props[clv_start:clv_end].mean(axis=0)
        return result
