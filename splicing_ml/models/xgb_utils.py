from __future__ import annotations

"""XGBoost GPU detection and inference-time device helpers.

Keeping GPU-related logic isolated here avoids scattering CUDA checks across
the pipeline and makes it easy to audit or disable GPU usage.
"""

import os
import subprocess
import warnings
from typing import Any

import numpy as np
from sklearn.pipeline import Pipeline

__all__ = [
    "xgb_gpu_available",
    "_set_xgb_cpu_predictor_for_inference",
    "_unwrap_model_step",
]

# Module-level cache so the CUDA probe runs at most once per process.
_XGB_GPU_AVAILABLE_CACHE: bool | None = None


def _unwrap_model_step(estimator: Any) -> Any:
    """Return the final model object when estimator is a sklearn Pipeline."""
    if isinstance(estimator, Pipeline) and "model" in estimator.named_steps:
        return estimator.named_steps["model"]
    return estimator


def _check_xgb_build_flags() -> bool | None:
    """Fast pre-check: return False if GPU is ruled out by build flags or env vars.

    Returns None when the check is inconclusive (runtime probe required).
    """
    try:
        import xgboost as xgb

        flag = str(xgb.build_info().get("USE_CUDA", "0")).lower()
        if flag not in {"1", "true", "on", "yes"}:
            return False
    except Exception:
        return False

    # Respect explicit CPU-forced runtime configuration.
    if str(os.environ.get("XGB_FORCE_CPU", "0")).lower() in {"1", "true", "yes", "on"}:
        return False

    # If the scheduler or container hides devices, treat as CPU-only.
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None:
        if visible.strip().lower() in {"", "-1", "none", "void", "n/a"}:
            return False

    return None  # inconclusive; runtime probe required


def _probe_cuda_runtime() -> bool:
    """Runtime CUDA probe: train a tiny XGBoost model on CUDA.

    A successful training run on device="cuda" is more reliable than
    environment-variable checks alone. Optionally cross-checks with
    nvidia-smi for diagnostic purposes (but never forces CPU-only based
    solely on that subprocess result).
    """
    try:
        import xgboost as xgb

        dprobe = xgb.DMatrix(
            np.asarray([[0.0], [1.0], [2.0], [3.0]], dtype=np.float32),
            label=np.asarray([0.0, 1.0, 0.0, 1.0], dtype=np.float32),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            booster = xgb.train(
                {
                    "objective": "binary:logistic",
                    "tree_method": "hist",
                    "device": "cuda",
                    "max_depth": 1,
                    "eta": 1.0,
                    "verbosity": 0,
                },
                dprobe,
                num_boost_round=1,
            )

        cfg_text = str(booster.save_config()).lower()
        if '"device":"cuda' not in cfg_text and '"device": "cuda' not in cfg_text:
            return False
    except Exception:
        return False

    # Best-effort diagnostics only; never force CPU solely because this
    # helper command is missing in containerised environments.
    try:
        subprocess.run(
            ["nvidia-smi", "-L"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=2,
            check=False,
        )
    except Exception:
        pass

    return True


def xgb_gpu_available() -> bool:
    """Best-effort CUDA runtime availability check for XGBoost.

    Requires both a CUDA-enabled XGBoost build and a visible runtime GPU.
    This avoids noisy warnings where XGBoost starts on CUDA and falls back
    because no GPU is actually available on the current node.

    Results are cached for the lifetime of the process.
    """
    global _XGB_GPU_AVAILABLE_CACHE
    if _XGB_GPU_AVAILABLE_CACHE is not None:
        return _XGB_GPU_AVAILABLE_CACHE

    fast_result = _check_xgb_build_flags()
    if fast_result is not None:
        _XGB_GPU_AVAILABLE_CACHE = fast_result
        return fast_result

    _XGB_GPU_AVAILABLE_CACHE = _probe_cuda_runtime()
    return _XGB_GPU_AVAILABLE_CACHE


def _set_xgb_cpu_predictor_for_inference(estimator: Any) -> None:
    """Switch fitted XGBoost estimators to CPU device before inference.

    The pipeline preprocessor outputs CPU arrays/matrices. If an XGBoost model
    remains on CUDA at predict-time, XGBoost emits repeated device-mismatch
    warnings and falls back internally. Setting ``device='cpu'`` after fitting
    avoids that mismatch while preserving GPU training in the tuning stage.
    """

    def _apply(obj: Any) -> None:
        if obj is None:
            return

        model = _unwrap_model_step(obj)
        model_module = getattr(model.__class__, "__module__", "")
        set_params = getattr(model, "set_params", None)

        # Only switch if model is actually on GPU to avoid unnecessary calls.
        if str(model_module).startswith("xgboost") and callable(set_params):
            try:
                current_device = getattr(model, "device", None)
                if current_device and str(current_device).lower() != "cpu":
                    set_params(device="cpu")
            except Exception:
                pass

        # Recurse into wrappers like CalibratedClassifierCV that hold children.
        for attr in ("estimator", "base_estimator", "classifier"):
            child = getattr(obj, attr, None)
            if child is not None and child is not obj:
                _apply(child)

        calibrated_items = getattr(obj, "calibrated_classifiers_", None)
        if calibrated_items is not None:
            for item in calibrated_items:
                _apply(item)

        estimators = getattr(obj, "estimators_", None)
        if estimators is not None:
            for item in estimators:
                _apply(item)

    _apply(estimator)
