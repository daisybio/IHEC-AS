from __future__ import annotations

"""GPU availability check and model factories for cuML estimators.

cuML (RAPIDS) provides GPU-accelerated estimators with near-identical APIs to
sklearn.  Detection is simpler than XGBoost: a bare ``import cuml`` suffices
because cuML triggers CUDA initialisation at import time, so an ImportError
definitively means no GPU support is available.

Supported cuML estimators:
- SVC / SVR                   (cuML kernel SVM, RBF by default; class_weight supported)
- RandomForestRegressor       (regression only; class_weight raises UnsupportedOnGPU)
- LinearRegression            (regression)
- LogisticRegression          (classification; supports class_weight="balanced")
- ElasticNet                  (regression; no CV-path variant in cuML)
- LogisticRegression          (classification elasticnet; penalty="elasticnet",
                               C and l1_ratio set via GridSearchCV externally)
"""

import os
from typing import Any

__all__ = [
    "cuml_gpu_available",
    "get_cuml_svm",
    "get_cuml_rf",
    "get_cuml_linear",
    "get_cuml_elasticnet",
    "get_cuml_logistic_elasticnet",
]

_CUML_GPU_AVAILABLE_CACHE: bool | None = None


def cuml_gpu_available() -> bool:
    """Return True when cuML is importable and a CUDA device is visible.

    Result is cached for the lifetime of the process (like ``xgb_gpu_available``).
    """
    global _CUML_GPU_AVAILABLE_CACHE
    if _CUML_GPU_AVAILABLE_CACHE is not None:
        return _CUML_GPU_AVAILABLE_CACHE

    # Respect explicit CPU override via environment variable.
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None and visible.strip().lower() in {
        "",
        "-1",
        "none",
        "void",
        "n/a",
    }:
        _CUML_GPU_AVAILABLE_CACHE = False
        return False

    try:
        import cuml  # noqa: F401

        _CUML_GPU_AVAILABLE_CACHE = True
    except Exception:
        _CUML_GPU_AVAILABLE_CACHE = False

    return _CUML_GPU_AVAILABLE_CACHE


def get_cuml_svm(task: str) -> Any:
    """Return a cuML SVC or SVR instance (RBF kernel by default).

    cuML SVC/SVR support C, gamma, and kernel parameters.  SVC also supports
    class_weight="balanced".  C and gamma are varied externally via GridSearchCV.
    """
    from cuml.svm import SVC, SVR

    if task == "classification":
        return SVC(kernel="rbf", class_weight="balanced", probability=True)
    return SVR(kernel="rbf")


def get_cuml_linear(task: str) -> Any:
    """Return a cuML LinearRegression or LogisticRegression instance.

    For classification, uses a very large C (1e6) to approximate the
    unregularized (C=inf) sklearn default.  class_weight="balanced" is
    supported by cuML LogisticRegression.
    """
    if task == "classification":
        from cuml.linear_model import LogisticRegression

        return LogisticRegression(C=1e6, class_weight="balanced", max_iter=5000)
    from cuml.linear_model import LinearRegression

    return LinearRegression()


def get_cuml_elasticnet() -> Any:
    """Return a cuML ElasticNet instance for regression.

    cuML ElasticNet accepts alpha and l1_ratio (same names as sklearn), but
    lacks the CV-path variant (ElasticNetCV).  The caller must handle
    hyperparameter search externally via GridSearchCV.
    """
    from cuml.linear_model import ElasticNet

    return ElasticNet(max_iter=20000)


def get_cuml_logistic_elasticnet() -> Any:
    """Return a base cuML LogisticRegression for elasticnet classification search.

    penalty="elasticnet" enables the combined L1+L2 penalty.  C and l1_ratio
    are left at defaults here and varied externally via GridSearchCV
    (model__C, model__l1_ratio).  class_weight="balanced" mirrors the sklearn
    LogisticRegressionCV fallback.
    """
    from cuml.linear_model import LogisticRegression

    return LogisticRegression(
        penalty="elasticnet",
        class_weight="balanced",
        max_iter=5000,
    )


def get_cuml_rf(task: str) -> Any:
    """Return a cuML RandomForestRegressor (regression only).

    cuML RandomForestClassifier raises UnsupportedOnGPU for class_weight,
    so classification falls back to the CPU sklearn path.

    cuML RF parameters are compatible with sklearn's RF, except:
      - max_features defaults to 1.0 instead of "auto"
    """
    if task == "classification":
        raise ValueError(
            "cuML RandomForestClassifier does not support class_weight parameter. "
            "Use CPU-based sklearn RandomForestClassifier instead."
        )
    from cuml.ensemble import RandomForestRegressor

    return RandomForestRegressor()
