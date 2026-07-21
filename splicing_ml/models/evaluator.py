from __future__ import annotations

"""Outer-fold evaluation helpers for regression and classification tasks.

Each task has a dedicated private function that handles its specific
logic (scale transformations, calibration, threshold tuning), with
``evaluate_outer_fold`` acting as a thin dispatcher.
"""

from typing import Any, Callable

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import Pipeline

from ..config import RNG_SEED
from ..metrics import (
    classification_metrics,
    regression_metrics,
    regression_metrics_by_psi_bin,
)
from ..utils import vlog
from .beta import BetaRegressor, _inverse_logit, _logit_transform
from .search import _lgbm_classification_sample_weight, _looks_like_gpu_failure
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


def _compute_xgb_shap(
    best_estimator: Pipeline,
    x_test: pd.DataFrame,
    verbose: bool,
    max_samples: int = 5000,
) -> dict[str, Any] | None:
    """Return per-test-sample SHAP values for an XGBoost/LightGBM outer-fold model.

    Detects XGBoost (sklearn API's ``get_booster`` method -- also present on
    XGBRFClassifier/Regressor, i.e. this project's "rf") and LightGBM (whose
    fitted sklearn wrapper exposes ``booster_`` instead, no ``get_booster``)
    — other model types return None. Uses ``best_estimator`` (the tuned
    pre-calibration Pipeline), never the Platt-calibrated wrapper:
    TreeExplainer needs the raw booster, not ``CalibratedClassifierCV``.

    max_samples: TreeExplainer's cost scales linearly with n_test_samples (no
    internal subsampling) -- on SE-scale folds (millions of rows) this took
    hours with zero progress output otherwise (found 2026-07-16 via a real
    ~2h13m silent gap in a production log). A random subsample of this size
    gives a stable mean-|SHAP| feature-importance estimate at a bounded,
    predictable cost regardless of fold size. ``sampled_test_indices`` in the
    return dict records which rows (positions into ``x_test``) were used, so
    results stay traceable back to the specific test rows. ``max_samples <= 0``
    disables subsampling entirely (use every test row) -- the intended
    setting for the final publication run, where the exact per-sample SHAP
    values matter more than wall-clock.
    """
    if not isinstance(best_estimator, Pipeline) or "model" not in best_estimator.named_steps:
        return None
    model = best_estimator.named_steps["model"]
    if not (hasattr(model, "get_booster") or hasattr(model, "booster_")):
        return None

    try:
        import time

        import shap

        n_test = len(x_test)
        if max_samples > 0 and n_test > max_samples:
            rng = np.random.default_rng(RNG_SEED)
            sample_pos = np.sort(rng.choice(n_test, size=max_samples, replace=False))
            x_shap = x_test.iloc[sample_pos]
        else:
            sample_pos = np.arange(n_test)
            x_shap = x_test

        prep = best_estimator.named_steps["prep"]
        x_test_t = prep.transform(x_shap)
        feature_names = (
            list(prep.get_feature_names_out())
            if hasattr(prep, "get_feature_names_out")
            else None
        )
        vlog(
            verbose,
            f"SHAP TreeExplainer starting: n_test={n_test}, n_sampled={len(x_test_t)}",
            level="info",
        )
        t0 = time.monotonic()
        explainer = shap.TreeExplainer(model)
        shap_values = np.asarray(explainer.shap_values(x_test_t))
        vlog(
            verbose,
            f"SHAP TreeExplainer done in {time.monotonic() - t0:.1f}s (n_sampled={len(x_test_t)})",
            level="info",
        )
        return {
            "shap_values": shap_values,
            "feature_names": feature_names,
            "sampled_test_indices": sample_pos,
        }
    except Exception as exc:
        compact = (
            str(exc).strip().splitlines()[0] if str(exc).strip() else exc.__class__.__name__
        )
        vlog(verbose, f"SHAP computation failed, skipping: {compact}", level="info")
        return None


def _compute_transformed_feature_sample(
    best_estimator: Pipeline,
    x_test: pd.DataFrame,
    verbose: bool,
    max_samples: int = 5000,
) -> dict[str, Any] | None:
    """Sample the fitted preprocessor's transformed output for distribution plots.

    Mirrors ``_compute_xgb_shap``'s subsampling (same ``RNG_SEED``, same
    bounded-cost rationale) but applies to any model -- this only needs the
    shared ``"prep"`` step's ``.transform()`` output, not a model-specific
    explainer. One-hot categorical columns (``"cat__"`` prefix, see
    ``build_preprocessor``) are excluded -- a 0/1 indicator's distribution
    isn't informative as a histogram.
    """
    if not isinstance(best_estimator, Pipeline) or "prep" not in best_estimator.named_steps:
        return None
    try:
        n_test = len(x_test)
        if max_samples > 0 and n_test > max_samples:
            rng = np.random.default_rng(RNG_SEED)
            sample_pos = np.sort(rng.choice(n_test, size=max_samples, replace=False))
            x_sample = x_test.iloc[sample_pos]
        else:
            x_sample = x_test

        prep = best_estimator.named_steps["prep"]
        x_t = prep.transform(x_sample)
        if not isinstance(x_t, pd.DataFrame):
            return None
        keep_cols = [c for c in x_t.columns if not str(c).startswith("cat__")]
        if not keep_cols:
            return None
        return {
            "feature_names": keep_cols,
            "values": x_t[keep_cols].to_numpy(dtype=float),
        }
    except Exception as exc:
        compact = (
            str(exc).strip().splitlines()[0] if str(exc).strip() else exc.__class__.__name__
        )
        vlog(
            verbose,
            f"Transformed-feature sampling failed, skipping: {compact}",
            level="info",
        )
        return None


def _chunked_model_call(
    call: Callable[[pd.DataFrame], np.ndarray],
    x: pd.DataFrame,
    verbose: bool,
    initial_chunk_size: int | None = None,
) -> np.ndarray:
    """Run a per-row model call (predict/predict_proba), shrinking into row-chunks only on OOM.

    TabICL is an in-context learner: predict_proba(X)/predict(X) take no
    batch-size argument, and its own internal batching (constructor-level
    ``batch_size``) doesn't adapt to memory pressure that appeared after
    fit() completed (e.g. a GPU that's gotten tighter since, or ordinary
    fragmentation) -- real SE-scale logs show it OOM-ing here even though
    fit() itself already succeeded (`Model failed: fold N model tabicl:
    CUDA out of memory` always follows a `Search complete model=tabicl`
    line).

    ``initial_chunk_size`` defaults to ``len(x)`` -- the first attempt is
    always a single call over the whole input, identical in cost to calling
    ``call(x)`` directly. Only an actual caught GPU OOM shrinks the chunk
    size and retries (down to a single row), with a cache clear between
    attempts. This matters because TabICL does not cache its fitted training
    context across separate predict calls -- each call re-materialises it
    from scratch (real logs show a ~90GB disk-offloaded rebuild per call) --
    so eagerly pre-chunking a call that would have succeeded whole multiplies
    that cost by the chunk count for no benefit (confirmed: a 522k-row test
    fold with a 50k starting chunk size took ~11x longer than a single call,
    turning a several-hour final eval into a ~16-hour one). Chunking is safe
    for every model type in this codebase (median-impute/scaler-based
    preprocessing and every model's predict/predict_proba are row-independent
    given a fixed fit), so this is applied unconditionally wherever it's
    wired in rather than gated to tabicl specifically -- the only cost for
    non-tabicl models on the happy path is now zero (single call, same as
    calling the model directly).
    """
    n = len(x)
    chunk_size = max(1, int(initial_chunk_size) if initial_chunk_size is not None else n)
    results: list[np.ndarray] = []
    start = 0
    while start < n:
        end = min(n, start + chunk_size)
        try:
            results.append(np.asarray(call(x.iloc[start:end])))
            start = end
        except Exception as exc:
            if not _looks_like_gpu_failure(exc) or chunk_size <= 1:
                raise
            chunk_size = max(1, chunk_size // 2)
            vlog(
                verbose,
                f"Chunked model call hit an apparent GPU OOM; shrinking "
                f"chunk_size to {chunk_size} and retrying ({exc})",
                level="info",
            )
            try:
                import torch

                torch.cuda.empty_cache()
            except ImportError:
                pass
    return np.concatenate(results, axis=0)


def tune_threshold_balanced_accuracy(
    estimator: Pipeline,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
) -> float:
    """Choose the decision threshold that maximises balanced accuracy.

    Predictions are collected over all inner validation folds and a grid of
    91 candidate thresholds [0.05, 0.95] (step 0.01) is evaluated. The
    threshold is locked here and must not be re-optimised on outer test data.

    Deliberately does NOT re-trigger early stopping on these per-inner-fold
    refits (fixed 2026-07-16): `estimator` arrives here already post-ES-refit
    (fit_best_estimator's tree_es step already fixed its real n_estimators,
    e.g. 89) -- `clone(estimator)` preserves that fixed count. Re-running
    early stopping here (as a prior version did, via _es_raw_val_context) let
    each of the ~4 per-fold refits independently stop at a *different* tree
    count than the model's real fixed n_estimators, on a further-shrunk
    training slice (fit on es_train_idx, an inner-fold cut of tr_idx) --
    producing out-of-fold probabilities from a family of inconsistently-sized
    models rather than a faithful proxy for the actual final model. Confirmed
    on real data (RI/High/seqnames, `.claude/scratch/
    lgbm_sample_weight_vs_scale_pos_weight_calibration.py`): this
    inconsistency was the dominant driver of a real, model-specific
    threshold-transfer gap between the chosen threshold and the outer-fold's
    actual balanced-accuracy-optimal threshold -- consistently worst for lgbm
    (leaf-wise growth, more sensitive to per-fold complexity swings) across
    every one of 5 real outer folds (avg gap 0.050, vs xgb's 0.020 and
    linear/tabicl's ~0.005-0.008, which have no analogous re-triggered-ES step
    at all). Fitting each inner-fold clone directly on its full `tr_idx` with
    the already-fixed n_estimators removes that inconsistency and uses
    slightly more data per fold besides.
    """
    thresholds = np.linspace(0.05, 0.95, 91)
    all_probs = np.zeros_like(y_train, dtype=float)
    # `estimator` here is the plain Pipeline fit_best_estimator returns after
    # _fit_tree_early_stopping's refit (search.best_estimator_ = pipe_es) --
    # NOT the _ESPipeline used during search/ES. XGBoost's scale_pos_weight is
    # baked into the constructor so it's applied automatically on every
    # .fit() call; LightGBM's imbalance correction was moved to a per-.fit()-
    # call sample_weight= (see _lgbm_classification_sample_weight) to dodge
    # the CUDA class-reweighting bug, so it does NOT survive clone()+plain
    # .fit() the way XGB's does -- it must be recomputed and passed here
    # explicitly, or every lgbm inner-fold refit silently trains with zero
    # imbalance correction. Found 2026-07-16 (second, distinct bug beyond the
    # ES-retriggering one fixed above): this was the DOMINANT remaining cause,
    # not inherent leaf-wise instability as first suspected -- real-data
    # verification across all 5 outer folds (RI/High/seqnames,
    # .claude/scratch/lgbm_threshold_fix_full_verify.py) moved the average
    # threshold-transfer gap 0.050 (original bug) -> 0.040 (ES-fix only) ->
    # 0.0063 (both fixes) -- on par with linear (0.005) and tabicl (0.008).
    lgbm_step = estimator.named_steps.get("model") if hasattr(estimator, "named_steps") else None
    is_lgbm = lgbm_step is not None and lgbm_step.__class__.__name__ == "LGBMClassifier"
    for tr_idx, va_idx in inner_splits:
        est = clone(estimator)
        if is_lgbm:
            sw = _lgbm_classification_sample_weight(y_train[tr_idx])
            est.fit(x_train.iloc[tr_idx], y_train[tr_idx], model__sample_weight=sw)
        else:
            est.fit(x_train.iloc[tr_idx], y_train[tr_idx])
        _set_xgb_cpu_predictor_for_inference(est)
        if hasattr(est, "predict_proba"):
            all_probs[va_idx] = _chunked_model_call(
                lambda xb: est.predict_proba(xb)[:, 1], x_train.iloc[va_idx], False
            )
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
    calibration_holdout: tuple[pd.DataFrame, np.ndarray] | None = None,
    shap_max_samples: int = 5000,
) -> dict[str, Any]:
    """Evaluate a classifier on one outer fold.

    Optionally tunes the decision threshold on inner validation predictions
    (never on outer test data), optionally calibrates probabilities with Platt
    scaling, and returns all classification metrics.

    calibration_holdout, when provided (XGB and LightGBM only — see
    ``_fit_tree_early_stopping``'s ``reserve_calibration_holdout``), is a
    (X, y) slice the fitted estimator has never seen; calibration then uses
    ``cv="prefit"`` on that slice instead of ``cv=inner_splits``, which would
    otherwise reuse the exact folds that already selected the hyperparameters.
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
        if calibration_holdout is not None:
            # Leakage-free path (XGB and LightGBM): calibrate on the ES-val slice the
            # fitted estimator never trained on, instead of cv=inner_splits
            # (which would reuse the folds that already picked the HPs).
            # sklearn >=1.6 removed cv="prefit"; FrozenEstimator is the
            # replacement for wrapping an already-fitted estimator.
            from sklearn.frozen import FrozenEstimator

            x_cal, y_cal = calibration_holdout
            calibrated = CalibratedClassifierCV(
                estimator=FrozenEstimator(best_estimator),
                method="sigmoid",
            )
            calibrated.fit(x_cal, y_cal)
        else:
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
            y_prob = _chunked_model_call(
                lambda xb: best_estimator.predict_proba(xb)[:, 1], x_test, verbose
            )
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
        "shap": _compute_xgb_shap(best_estimator, x_test, verbose, shap_max_samples),
        "transformed_features": _compute_transformed_feature_sample(
            best_estimator, x_test, verbose, shap_max_samples
        ),
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
    shap_max_samples: int = 5000,
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
    y_pred_model = _chunked_model_call(best_estimator.predict, x_test, verbose).astype(
        float
    )

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
        "shap": _compute_xgb_shap(best_estimator, x_test, verbose, shap_max_samples),
        "transformed_features": _compute_transformed_feature_sample(
            best_estimator, x_test, verbose, shap_max_samples
        ),
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
    calibration_holdout: tuple[pd.DataFrame, np.ndarray] | None = None,
    shap_max_samples: int = 5000,
) -> dict[str, Any]:
    """Evaluate a tuned estimator on one outer fold.

    Dispatches to ``_eval_classification`` or ``_eval_regression`` based on
    task. Both return a dict with keys: scores, threshold, y_true, y_pred,
    estimator, shap (plus y_true_logit / y_pred_logit for regression). ``shap``
    is a ``{"shap_values", "feature_names", "sampled_test_indices"}`` dict for
    XGB/LightGBM models, else None.

    Parameters
    ----------
    calibrate
        Whether to apply Platt probability calibration for classifiers.
        Calibration is slower but generally recommended.
    tune_threshold
        Whether to search for the optimal decision threshold on inner-fold
        predictions.  When False, threshold is fixed at 0.5.
    calibration_holdout
        Optional (X, y) slice unseen by ``best_estimator`` (XGB and LightGBM
        only, from ``tuning_info["calibration_holdout"]``); when present,
        calibration wraps ``best_estimator`` in ``FrozenEstimator`` and fits
        on it instead of using ``cv=inner_splits``.
    shap_max_samples
        Max outer-test rows used for SHAP TreeExplainer (xgb/lgbm only) --
        see ``_compute_xgb_shap``.
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
            calibration_holdout=calibration_holdout,
            shap_max_samples=shap_max_samples,
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
        shap_max_samples=shap_max_samples,
    )
