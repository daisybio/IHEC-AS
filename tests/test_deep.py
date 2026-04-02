from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse
from sklearn.preprocessing import StandardScaler

from splicing_ml.models._contexts import _es_pp_val_context
from splicing_ml.models.deep import MLPRegressor
from splicing_ml.models.search import _ESPipeline, _es_raw_val_context


def test_mlp_regressor_handles_sparse_validation_context() -> None:
    """Regression test for sparse preprocessed validation data."""
    torch = pytest.importorskip("torch")

    rng = np.random.default_rng(0)
    x_train = rng.normal(size=(8, 4)).astype(np.float32)
    y_train = rng.normal(size=8).astype(np.float32)
    x_val = sparse.csr_matrix(rng.normal(size=(3, 4)).astype(np.float32))
    y_val = rng.normal(size=3).astype(np.float32)

    _es_pp_val_context.X_val_pp = x_val
    _es_pp_val_context.y_val = y_val
    try:
        model = MLPRegressor(
            hidden_sizes=(4,),
            dropout_rate=0.0,
            learning_rate=1e-3,
            weight_decay=0.0,
            batch_size=2,
            max_epochs=1,
            early_stopping_patience=1,
            random_state=0,
        )
        fitted = model.fit(x_train, y_train)

        assert fitted is model
        assert getattr(model, "n_epochs_trained_", 0) >= 1
    finally:
        _es_pp_val_context.X_val_pp = None
        _es_pp_val_context.y_val = None


@pytest.mark.parametrize("model_name", ["mlp", "xgb"])
def test_es_pipeline_handles_validation_routing(model_name: str) -> None:
    """Regression test for sklearn metadata routing in the early-stopping pipeline."""
    torch = pytest.importorskip("torch")

    rng = np.random.default_rng(1)
    x_train = rng.normal(size=(6, 4)).astype(np.float32)
    y_train = rng.normal(size=6).astype(np.float32)
    x_val = rng.normal(size=(2, 4)).astype(np.float32)
    y_val = rng.normal(size=2).astype(np.float32)

    if model_name == "mlp":
        model = MLPRegressor(
            hidden_sizes=(4,),
            dropout_rate=0.0,
            learning_rate=1e-3,
            weight_decay=0.0,
            batch_size=2,
            max_epochs=1,
            early_stopping_patience=1,
            random_state=0,
        )
    else:
        xgboost = pytest.importorskip("xgboost")
        model = xgboost.XGBRegressor(
            n_estimators=2,
            max_depth=2,
            learning_rate=0.1,
            subsample=1.0,
            colsample_bytree=1.0,
            min_child_weight=1,
            tree_method="hist",
            device="cpu",
            eval_metric="rmse",
            random_state=0,
        )

    pipeline = _ESPipeline(
        steps=[("prep", StandardScaler()), ("model", model)],
        model_name=model_name,
    )

    _es_raw_val_context.X_val = x_val
    _es_raw_val_context.y_val = y_val
    try:
        fitted = pipeline.fit(x_train, y_train)
        assert fitted is pipeline
    finally:
        _es_raw_val_context.X_val = None
        _es_raw_val_context.y_val = None
