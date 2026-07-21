from __future__ import annotations

"""LightGBM GPU detection and inference-time device helpers.

Mirrors xgb_utils.py's shape for the boosted-tree model family, but LightGBM
is a separate library from XGBoost with its own GPU mechanism and no
`build_info()`-style API, so the probe strategy differs (env-flag pre-check +
a real runtime fit, no build-flag short-circuit).
"""

import os
import warnings
from typing import Any

import numpy as np

from .xgb_utils import _unwrap_model_step  # fully generic, reused as-is

__all__ = [
    "lgbm_gpu_available",
    "_set_lgbm_cpu_predictor_for_inference",
]

# Module-level cache so the CUDA probe runs at most once per process.
_LGBM_GPU_AVAILABLE_CACHE: bool | None = None


def _check_lgbm_env_flags() -> bool | None:
    """Fast pre-check: return False if GPU is ruled out by env vars.

    Unlike XGBoost, LightGBM has no `build_info()`-style API to check whether
    the installed build even supports CUDA, so this only handles the
    explicit-override and container-hides-devices cases; everything else
    falls through to the runtime probe.
    """
    if str(os.environ.get("LGBM_FORCE_CPU", "0")).lower() in {"1", "true", "yes", "on"}:
        return False

    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None and visible.strip().lower() in {"", "-1", "none", "void", "n/a"}:
        return False

    return None  # inconclusive; runtime probe required


def _probe_lgbm_cuda_runtime() -> bool:
    """Runtime CUDA probe: fit device_type='cuda' on separable synthetic data,
    using the same eval_set + early_stopping callback shape production fits
    always use, and check it actually learns (not just that .fit() runs).

    A crash-only probe (the original version of this function) is not
    sufficient: a real run on 2026-07-14 hit device_type="cuda" fitting
    without error yet producing exact-chance predictions (AUROC=0.500 on
    every one of 16 hyperparameter candidates) on real data, while an
    otherwise-identical device_type="cpu" fit reached AUROC=0.58-0.70 on the
    same data/features through the exact same _ESPipeline/grid-search code
    path — confirmed via a real-pipeline CPU repro, so this is CPU-vs-GPU
    device-specific, not a plumbing bug. A first version of this probe (fit
    with no eval_set/callback at all) still reported GPU as healthy on real
    hardware where the actual failure reproduced, because every real lgbm
    fit in this codebase always passes eval_set + a
    lightgbm.early_stopping(...) callback (see _ESPipeline.fit and
    _fit_tree_early_stopping) — device_type="cuda" combined with
    callback-based early stopping is a known trouble spot for LightGBM's CUDA
    backend, and a probe that skips that combination can't catch it. So this
    probe now mirrors production exactly: eval_set + early_stopping callback,
    device_type="cuda", and requires a near-perfect held-out AUROC.
    """
    try:
        import lightgbm
        from lightgbm import LGBMClassifier
        from sklearn.metrics import roc_auc_score

        rng = np.random.default_rng(0)
        x = rng.normal(size=(200, 4)).astype(np.float32)
        y = (x[:, 0] + 0.5 * x[:, 1] > 0).astype(int)
        x_train, y_train = x[:120], y[:120]
        x_val, y_val = x[120:160], y[120:160]
        x_test, y_test = x[160:], y[160:]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = LGBMClassifier(
                device_type="cuda",
                n_estimators=50,
                num_leaves=7,
                min_child_samples=1,
                verbose=-1,
                random_state=0,
            )
            model.fit(
                x_train,
                y_train,
                eval_set=[(x_val, y_val)],
                callbacks=[lightgbm.early_stopping(stopping_rounds=10, verbose=False)],
            )
            preds = model.predict_proba(x_test)[:, 1]

        if len(set(y_test.tolist())) < 2:
            return False
        auroc = roc_auc_score(y_test, preds)
    except Exception:
        return False

    return auroc > 0.9


def lgbm_gpu_available() -> bool:
    """Best-effort CUDA runtime *and correctness* check for LightGBM.

    Only probes the dedicated CUDA device type ("cuda"), not the OpenCL one
    ("gpu") — matches the cuda129-tagged conda-forge build this project pins.
    Beyond "does it crash", also checks the fit actually learns (see
    `_probe_lgbm_cuda_runtime`'s docstring) — a crash-only check previously
    let a silently-broken cuda backend through. Results are cached for the
    lifetime of the process.
    """
    global _LGBM_GPU_AVAILABLE_CACHE
    if _LGBM_GPU_AVAILABLE_CACHE is not None:
        return _LGBM_GPU_AVAILABLE_CACHE

    fast_result = _check_lgbm_env_flags()
    if fast_result is not None:
        _LGBM_GPU_AVAILABLE_CACHE = fast_result
        return fast_result

    _LGBM_GPU_AVAILABLE_CACHE = _probe_lgbm_cuda_runtime()
    return _LGBM_GPU_AVAILABLE_CACHE


def _set_lgbm_cpu_predictor_for_inference(estimator: Any) -> None:
    """Switch fitted LightGBM estimators to CPU device before inference.

    Structurally identical to xgb_utils._set_xgb_cpu_predictor_for_inference,
    but checks the "lightgbm" module prefix and LightGBM's `device_type`
    constructor param (not XGBoost's `device`).
    """

    def _apply(obj: Any) -> None:
        """Internal helper for apply."""
        if obj is None:
            return

        model = _unwrap_model_step(obj)
        model_module = getattr(model.__class__, "__module__", "")
        set_params = getattr(model, "set_params", None)

        if str(model_module).startswith("lightgbm") and callable(set_params):
            try:
                current_device = getattr(model, "device_type", None)
                if current_device and str(current_device).lower() != "cpu":
                    set_params(device_type="cpu")
            except Exception:
                pass

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
