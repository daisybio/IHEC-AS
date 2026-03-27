from __future__ import annotations

"""GPU availability check and model factory for cuML SVM.

cuML (RAPIDS) provides GPU-accelerated LinearSVC/LinearSVR with near-identical
APIs to sklearn.  Detection is simpler than XGBoost: a bare ``import cuml``
suffices because cuML triggers CUDA initialisation at import time, so an
ImportError definitively means no GPU support is available.

cuML uses an ADMM GPU solver and does not accept a ``dual`` parameter.
"""

import os
from typing import Any

__all__ = ["cuml_gpu_available", "get_cuml_svm"]

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
    if visible is not None and visible.strip().lower() in {"", "-1", "none", "void", "n/a"}:
        _CUML_GPU_AVAILABLE_CACHE = False
        return False

    try:
        import cuml  # noqa: F401

        _CUML_GPU_AVAILABLE_CACHE = True
    except Exception:
        _CUML_GPU_AVAILABLE_CACHE = False

    return _CUML_GPU_AVAILABLE_CACHE


def get_cuml_svm(task: str) -> Any:
    """Return a cuML LinearSVC or LinearSVR instance with best-practice defaults.

    cuML's ADMM GPU solver does not accept a ``dual`` parameter (unlike sklearn's
    liblinear-backed LinearSVC/LinearSVR).  All other parameters match the
    sklearn fallback used in grids.py.
    """
    from cuml.svm import LinearSVC, LinearSVR

    if task == "classification":
        return LinearSVC(class_weight="balanced", max_iter=5000)
    return LinearSVR(epsilon=0.0, max_iter=5000)
