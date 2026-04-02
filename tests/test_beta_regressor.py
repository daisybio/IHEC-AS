from __future__ import annotations

import numpy as np

from splicing_ml.models.beta import BetaRegressor


class _DummyBetaModel:
    def __init__(self, n_cols: int) -> None:
        """Initialize a _DummyBetaModel instance."""
        self.exog = np.zeros((1, n_cols), dtype=float)


class _DummyBetaResult:
    def __init__(self, n_cols: int) -> None:
        """Initialize a _DummyBetaResult instance."""
        self.model = _DummyBetaModel(n_cols)

    def predict(self, x_design: np.ndarray) -> np.ndarray:
        """Generate predictions for the provided input data."""
        if x_design.shape[1] != self.model.exog.shape[1]:
            raise ValueError(
                f"shape mismatch: got {x_design.shape[1]}, expected {self.model.exog.shape[1]}"
            )
        return np.full(x_design.shape[0], 0.5, dtype=float)


def test_predict_adds_intercept_even_with_constant_feature() -> None:
    """Regression test for fold-dependent intercept drops in predict()."""
    reg = BetaRegressor()
    # Fitted model expects 49 columns total (48 features + intercept).
    reg._beta_result_ = _DummyBetaResult(n_cols=49)

    x = np.random.RandomState(0).randn(8, 48)
    # Constant feature can trigger add_constant(..., has_constant='skip') behavior.
    x[:, 0] = 1.0

    y_pred = reg.predict(x)
    assert y_pred.shape == (8,)
    assert np.all((y_pred > 0.0) & (y_pred < 1.0))
