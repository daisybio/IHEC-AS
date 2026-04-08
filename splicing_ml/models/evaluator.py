from __future__ import annotations

"""Outer-fold evaluation helpers for regression and classification tasks.

Each task has a dedicated private function that handles its specific
logic (scale transformations, calibration, threshold tuning), with
``evaluate_outer_fold`` acting as a thin dispatcher.
"""

from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import Pipeline

from ..metrics import (
    classification_metrics,
    regression_metrics,
    regression_metrics_by_psi_bin,
)
from ..utils import vlog
from ._contexts import _es_raw_val_context
from .beta import BetaRegressor, _inverse_logit, _logit_transform
from .search import _compute_es_splits
from .xgb_utils import _set_xgb_cpu_predictor_for_inference

__all__ = ["evaluate_outer_fold", "tune_threshold_balanced_accuracy"]


def _model_outputs_psi_scale(estimator: Any) -> bool:
    """Return True when the estimator predicts directly on the PSI [0, 1] scale."""
    model = (
        estimator.named_steps["model"]
        if isinstance(estimator, Pipeline) and "model" in estimator.named_steps
        else estimator
    )
    return isinstance(model, BetaRegressor)


def tune_threshold_balanced_accuracy(
    estimator: Pipeline,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
) -> float:
    """Choose the decision threshold that maximises balanced accuracy.

    Predictions are collected over all inner validation folds and a grid of
    19 candidate thresholds [0.05, 0.95] is evaluated. The threshold is
    locked here and must not be re-optimised on outer test data.
    """
    thresholds = np.linspace(0.05, 0.95, 19)
    all_probs = np.zeros_like(y_train, dtype=float)
    # Mirror _ContextAwareCV/_compute_es_splits: dedicated ES val from another fold,
    # disjoint from the scoring val. Only computed when >=2 splits available.
    es_splits = _compute_es_splits(inner_splits) if len(inner_splits) >= 2 else None
    for i, (tr_idx, va_idx) in enumerate(inner_splits):
        est = clone(estimator)
        if es_splits is not None:
            es_train_idx, es_val_idx = es_splits[i]
            _es_raw_val_context.X_val = x_train.iloc[es_val_idx]
            _es_raw_val_context.y_val = y_train[es_val_idx]
            fit_idx = es_train_idx
        else:
            fit_idx = tr_idx
        try:
            est.fit(x_train.iloc[fit_idx], y_train[fit_idx])
        finally:
            _es_raw_val_context.X_val = None
            _es_raw_val_context.y_val = None
        _set_xgb_cpu_predictor_for_inference(est)
        if hasattr(est, "predict_proba"):
            all_probs[va_idx] = est.predict_proba(x_train.iloc[va_idx])[:, 1]
        else:
            decision = est.decision_function(x_train.iloc[va_idx])
            all_probs[va_idx] = 1.0 / (1.0 + np.exp(-decision))

    best_threshold, best_score = 0.5, -np.inf
    for thr in thresholds:
        score = balanced_accuracy_score(y_train, (all_probs >= thr).astype(int))
        if score > best_score:
            best_score, best_threshold = score, float(thr)
    return best_threshold


def _eval_classification(
    best_estimator: Pipeline,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    x_test: pd.DataFrame,
    y_test: np.ndarray,
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    calibrate: bool,
    tune_threshold: bool,
    verbose: bool,
) -> dict[str, Any]:
    """Evaluate a classifier on one outer fold.

    Optionally tunes the decision threshold on inner validation predictions
    (never on outer test data), optionally calibrates probabilities with Platt
    scaling, and returns all classification metrics.
    """
    if tune_threshold:
        threshold = tune_threshold_balanced_accuracy(
            best_estimator, x_train, y_train, inner_splits
        )
        vlog(verbose, f"Threshold tuning selected threshold={threshold:.4f}")
    else:
        threshold = 0.5
        vlog(verbose, f"Threshold tuning disabled; threshold={threshold:.4f}")

    if calibrate:
        # Fit Platt-scaled calibration using the inner splits as CV folds.
        calibrated = CalibratedClassifierCV(
            estimator=best_estimator,
            method="sigmoid",
            cv=inner_splits,
            n_jobs=max_cores,
        )
        calibrated.fit(x_train, y_train)
        _set_xgb_cpu_predictor_for_inference(calibrated)
        y_prob = calibrated.predict_proba(x_test)[:, 1]
        final_estimator = calibrated
    else:
        _set_xgb_cpu_predictor_for_inference(best_estimator)
        if hasattr(best_estimator, "predict_proba"):
            y_prob = best_estimator.predict_proba(x_test)[:, 1]
        else:
            decision = best_estimator.decision_function(x_test)
            y_prob = 1.0 / (1.0 + np.exp(-decision))
        final_estimator = best_estimator

    scores = classification_metrics(y_test, y_prob, threshold=threshold)
    return {
        "scores": scores,
        "threshold": threshold,
        "y_true": y_test,
        "y_pred": y_prob,
        "estimator": final_estimator,
    }


def _eval_regression(
    best_estimator: Pipeline,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    x_test: pd.DataFrame,
    y_test: np.ndarray,
    verbose: bool,
    psi_low: float = 0.2,
    psi_high: float = 0.8,
) -> dict[str, Any]:
    """Evaluate a regressor on one outer fold.

    Handles the dual-scale reporting requirement: metrics are computed on
    both the original PSI scale and the logit scale, and both sets are
    exposed in the return dict (prefixed ``original_`` and ``logit_``).
    The legacy (un-prefixed) metric names map to the original PSI scale
    for backward compatibility.
    """
    y_train_fit = np.asarray(y_train, dtype=float)
    y_test_eval = np.asarray(y_test, dtype=float)

    # BetaRegressor expects raw PSI targets; invert logit if needed.
    if _model_outputs_psi_scale(best_estimator) and (
        np.any(y_train_fit < 0.0) or np.any(y_train_fit > 1.0)
    ):
        y_train_fit = _inverse_logit(y_train_fit)

    best_estimator.fit(x_train, y_train_fit)
    _set_xgb_cpu_predictor_for_inference(best_estimator)
    y_pred_model = np.asarray(best_estimator.predict(x_test), dtype=float)

    # Determine whether the test target is on the logit or PSI scale and
    # align predictions accordingly before computing metrics.
    is_logit_target = bool(np.any(y_test_eval < 0.0) or np.any(y_test_eval > 1.0))
    if is_logit_target:
        y_test_original = _inverse_logit(y_test_eval)
        y_pred_original = (
            np.clip(y_pred_model, 0.0, 1.0)
            if _model_outputs_psi_scale(best_estimator)
            else _inverse_logit(y_pred_model)
        )
    else:
        y_test_original = y_test_eval
        y_pred_original = y_pred_model

    scores_original = regression_metrics(y_test_original, y_pred_original)
    y_test_logit = _logit_transform(y_test_original)
    y_pred_logit = _logit_transform(y_pred_original)
    scores_logit = regression_metrics(y_test_logit, y_pred_logit)

    # Build combined scores dict: legacy un-prefixed names (PSI scale) +
    # explicit original_/logit_ prefixed variants for dual-scale reporting.
    scores = dict(scores_original)
    for k, v in scores_original.items():
        scores[f"original_{k}"] = float(v)
    for k, v in scores_logit.items():
        scores[f"logit_{k}"] = float(v)

    psi_bin_metrics = regression_metrics_by_psi_bin(
        y_test_original, y_pred_original, psi_low=psi_low, psi_high=psi_high
    )

    return {
        "scores": scores,
        "threshold": None,
        "y_true": y_test_original,
        "y_pred": y_pred_original,
        "y_true_logit": y_test_logit,
        "y_pred_logit": y_pred_logit,
        "estimator": best_estimator,
        "psi_bin_metrics": psi_bin_metrics,
    }


def evaluate_outer_fold(
    task: str,
    best_estimator: Pipeline,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    x_test: pd.DataFrame,
    y_test: np.ndarray,
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    calibrate: bool = True,
    tune_threshold: bool = True,
    verbose: bool = False,
    psi_low: float = 0.2,
    psi_high: float = 0.8,
) -> dict[str, Any]:
    """Evaluate a tuned estimator on one outer fold.

    Dispatches to ``_eval_classification`` or ``_eval_regression`` based on
    task. Both return a dict with keys: scores, threshold, y_true, y_pred,
    estimator (plus y_true_logit / y_pred_logit for regression).

    Parameters
    ----------
    calibrate
        Whether to apply Platt probability calibration for classifiers.
        Calibration is slower but generally recommended.
    tune_threshold
        Whether to search for the optimal decision threshold on inner-fold
        predictions.  When False, threshold is fixed at 0.5.
    """
    if task == "classification":
        return _eval_classification(
            best_estimator=best_estimator,
            x_train=x_train,
            y_train=y_train,
            x_test=x_test,
            y_test=y_test,
            inner_splits=inner_splits,
            max_cores=max_cores,
            calibrate=calibrate,
            tune_threshold=tune_threshold,
            verbose=verbose,
        )
    return _eval_regression(
        best_estimator=best_estimator,
        x_train=x_train,
        y_train=y_train,
        x_test=x_test,
        y_test=y_test,
        verbose=verbose,
        psi_low=psi_low,
        psi_high=psi_high,
    )
