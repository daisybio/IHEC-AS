from __future__ import annotations

"""Plot payload builders for HTML reports.

The main entry point is ``build_task_plot_payload``, which converts a single
task result dict into a JSON-serialisable plotting payload consumed by the
Plotly-based JavaScript renderer in the HTML report.

Key improvements over the original monolithic function:
- ``_build_metric_heatmap`` promoted from nested to module-level (called 3×).
- ``_safe_empty_payload`` extracts the two early-return edge cases.
- ``_build_predictions_payload`` and ``_build_classification_threshold_payload``
  separate the main prediction/threshold/confusion collection loop.
"""

import collections
import json
from typing import Any

import numpy as np

from ..metrics import regression_metrics
from ..utils import safe_json

__all__ = ["build_task_plot_payload", "important_params_table_rows"]


def _safe_logit(values: list[float], epsilon: float = 1e-7) -> list[float]:
    """Convert PSI-like values to logit scale with epsilon clipping."""
    if not values:
        return []
    arr = np.asarray(values, dtype=float)
    arr = np.clip(arr, epsilon, 1.0 - epsilon)
    return np.log(arr / (1.0 - arr)).tolist()


def important_params_table_rows(task_result: dict[str, Any]) -> list[dict[str, str]]:
    """Summarise the most frequently selected best-parameter combinations per model."""
    by_model: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in task_result.get("fold_results", []):
        by_model[str(row.get("model_name", "unknown"))].append(
            row.get("tuning", {}).get("best_params", {})
        )

    out: list[dict[str, str]] = []
    for model_name, params in by_model.items():
        ctr = collections.Counter(
            json.dumps(safe_json(p), sort_keys=True) for p in params
        )
        if not ctr:
            continue
        p_json, freq = ctr.most_common(1)[0]
        out.append(
            {
                "model_name": model_name,
                "selected_in_folds": str(len(params)),
                "most_frequent_best_params": p_json,
                "frequency": str(freq),
            }
        )
    return out


# ---------------------------------------------------------------------------
# Heatmap builder (promoted from nested function — called 3× in payload)
# ---------------------------------------------------------------------------


def _build_metric_heatmap(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
    metric_names: list[str],
) -> dict[str, Any]:
    """Build a Plotly-ready heatmap dict with means, SDs, and formatted text.

    Parameters
    ----------
    fold_rows : list of fold result dicts
    model_order : model names in display order (y-axis)
    metric_names : metric keys to include (x-axis)
    """
    metric_z: list[list[float | None]] = []
    metric_sd: list[list[float | None]] = []
    metric_text: list[list[str]] = []
    for model_name in model_order:
        model_rows = [
            r for r in fold_rows if str(r.get("model_name", "unknown")) == model_name
        ]
        row_means: list[float | None] = []
        row_sds: list[float | None] = []
        row_text: list[str] = []
        for metric in metric_names:
            vals = [
                float(r.get("scores", {}).get(metric, np.nan))
                for r in model_rows
                if not np.isnan(float(r.get("scores", {}).get(metric, np.nan)))
            ]
            if vals:
                mean_val = float(np.mean(vals))
                sd_val = float(np.std(vals))
                row_means.append(mean_val)
                row_sds.append(sd_val)
                row_text.append(f"{mean_val:.3f} +/- {sd_val:.3f}")
            else:
                row_means.append(None)
                row_sds.append(None)
                row_text.append("NA")
        metric_z.append(row_means)
        metric_sd.append(row_sds)
        metric_text.append(row_text)
    return {
        "models": model_order,
        "metrics": metric_names,
        "z": metric_z,
        "sd": metric_sd,
        "text": metric_text,
    }


_EMPTY_HEATMAP: dict[str, Any] = {
    "models": [],
    "metrics": [],
    "z": [],
    "sd": [],
    "text": [],
}


def _safe_empty_payload(
    task_result: dict[str, Any],
    model_fold_points: list[dict[str, Any]] | None,
    reason: str,
) -> dict[str, Any]:
    """Return a minimal safe payload when fold results are absent or all-NaN."""
    warnings = list(task_result.get("warnings", []))
    warnings.append(reason)
    return {
        "task": task_result.get("task", "unknown"),
        "status": task_result.get("status", "unknown"),
        "primary_metric": task_result.get("primary_metric", "metric"),
        "n_samples": task_result.get("n_samples", None),
        "warnings": warnings,
        "model_order": [],
        "metric_heatmap": _EMPTY_HEATMAP.copy(),
        "model_fold_points": model_fold_points or [],
        "fold_lines": [],
        "y_true": [],
        "y_pred": [],
        "prediction_by_model": {},
        "prediction_by_model_scale": {},
        "thresholds": [],
        "threshold_by_model": [],
        "best_threshold": None,
        "confusion_matrix": None,
        "confusion_by_model": {},
        "roc_by_model": {},
        "pr_by_model": {},
        "response_distribution": _ensure_thresholds(task_result),
        "important_params": important_params_table_rows(task_result),
    }


# ---------------------------------------------------------------------------
# Response-distribution helpers
# ---------------------------------------------------------------------------


def _ensure_thresholds(task_result: dict[str, Any]) -> dict[str, Any]:
    """Return response_distribution with binarization_thresholds always populated.

    Older result JSONs store ``null`` for regression. Fall back to
    ``reproducibility.psi_thresholds`` (always written by the pipeline) so
    the HTML report can shade excluded regions for both tasks.
    """
    dist = dict(task_result.get("response_distribution", {}))
    if not dist.get("binarization_thresholds"):
        fallback = (task_result.get("reproducibility") or {}).get("psi_thresholds")
        if fallback:
            dist["binarization_thresholds"] = fallback
    return dist


# ---------------------------------------------------------------------------
# Prediction and threshold payload builders
# ---------------------------------------------------------------------------


def _build_predictions_payload(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
    task: str,
) -> tuple[
    list[float],
    list[float],
    list[float],
    dict[str, dict[str, list[float]]],
    dict[str, dict[str, dict[str, list[float]]]],
]:
    """Collect y_true, y_pred, thresholds, and per-model prediction dicts."""
    y_true: list[float] = []
    y_pred: list[float] = []
    thresholds: list[float] = []
    prediction_by_model: dict[str, dict[str, list[float]]] = {
        m: {"y_true": [], "y_pred": [], "thresholds": []} for m in model_order
    }
    prediction_by_model_scale: dict[str, dict[str, dict[str, list[float]]]] = {
        m: {
            "original": {"y_true": [], "y_pred": []},
            "logit": {"y_true": [], "y_pred": []},
        }
        for m in model_order
    }

    for r in fold_rows:
        row_model = str(r.get("model_name", "unknown"))
        row_y_true = r.get("y_true", [])
        row_y_pred = r.get("y_pred", [])
        row_y_true_logit = r.get("y_true_logit", [])
        row_y_pred_logit = r.get("y_pred_logit", [])

        # Compute logit-scale values on the fly if not stored (older artifacts).
        if task == "regression":
            if not row_y_true_logit and row_y_true:
                row_y_true_logit = _safe_logit(row_y_true)
            if not row_y_pred_logit and row_y_pred:
                row_y_pred_logit = _safe_logit(row_y_pred)

        y_true.extend(row_y_true)
        y_pred.extend(row_y_pred)

        thr = r.get("threshold")
        thr_val = float(thr) if thr is not None else 0.5
        thresholds.append(thr_val)

        if row_model in prediction_by_model:
            prediction_by_model[row_model]["y_true"].extend(row_y_true)
            prediction_by_model[row_model]["y_pred"].extend(row_y_pred)
            prediction_by_model[row_model]["thresholds"].append(thr_val)

        if row_model in prediction_by_model_scale:
            prediction_by_model_scale[row_model]["original"]["y_true"].extend(row_y_true)
            prediction_by_model_scale[row_model]["original"]["y_pred"].extend(row_y_pred)
            prediction_by_model_scale[row_model]["logit"]["y_true"].extend(row_y_true_logit)
            prediction_by_model_scale[row_model]["logit"]["y_pred"].extend(row_y_pred_logit)

    # Down-sample to keep HTML lightweight for large runs.
    if len(y_true) > 8000:
        idx = np.linspace(0, len(y_true) - 1, num=8000, dtype=int)
        y_true = [y_true[i] for i in idx]
        y_pred = [y_pred[i] for i in idx]

    return y_true, y_pred, thresholds, prediction_by_model, prediction_by_model_scale


def _build_roc_payload(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
) -> dict[str, dict[str, Any]]:
    """Build per-model ROC curve data by pooling predictions across folds.

    Returns a dict keyed by model name, each value containing
    ``fpr``, ``tpr`` (lists), and ``auc`` (float).  Models with
    fewer than two classes in the pooled labels are skipped.
    """
    from sklearn.metrics import auc, roc_curve

    roc_by_model: dict[str, dict[str, Any]] = {}
    for model_name in model_order:
        m_rows = [r for r in fold_rows if str(r.get("model_name", "")) == model_name]
        all_y_true: list[float] = []
        all_y_pred: list[float] = []
        for r in m_rows:
            all_y_true.extend(r.get("y_true", []))
            all_y_pred.extend(r.get("y_pred", []))
        if not all_y_true:
            continue
        yt = np.asarray(all_y_true, dtype=int)
        yp = np.asarray(all_y_pred, dtype=float)
        if len(np.unique(yt)) < 2:
            continue
        fpr, tpr, _ = roc_curve(yt, yp)
        auc_val = float(auc(fpr, tpr))
        roc_by_model[model_name] = {
            "fpr": fpr.tolist(),
            "tpr": tpr.tolist(),
            "auc": auc_val,
        }
    return roc_by_model


def _build_pr_payload(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
) -> dict[str, dict[str, Any]]:
    """Build per-model PR curve data by pooling predictions across folds.

    Returns a dict keyed by model name, each value containing
    ``precision``, ``recall`` (lists), and ``avg_precision`` (float).
    Models with fewer than two classes in the pooled labels are skipped.
    """
    from sklearn.metrics import average_precision_score, precision_recall_curve

    pr_by_model: dict[str, dict[str, Any]] = {}
    for model_name in model_order:
        m_rows = [r for r in fold_rows if str(r.get("model_name", "")) == model_name]
        all_y_true: list[float] = []
        all_y_pred: list[float] = []
        for r in m_rows:
            all_y_true.extend(r.get("y_true", []))
            all_y_pred.extend(r.get("y_pred", []))
        if not all_y_true:
            continue
        yt = np.asarray(all_y_true, dtype=int)
        yp = np.asarray(all_y_pred, dtype=float)
        if len(np.unique(yt)) < 2:
            continue
        precision, recall, _ = precision_recall_curve(yt, yp)
        ap = float(average_precision_score(yt, yp))
        pr_by_model[model_name] = {
            "precision": precision.tolist(),
            "recall": recall.tolist(),
            "avg_precision": ap,
        }
    return pr_by_model


def _build_classification_threshold_payload(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
    primary_metric: str,
) -> tuple[
    dict[str, list[list[int]]],
    list[dict[str, Any]],
    dict[str, Any] | None,
]:
    """Compute per-model confusion matrices and threshold summary statistics."""
    confusion_by_model: dict[str, list[list[int]]] = {}
    threshold_summary: list[dict[str, Any]] = []

    for model_name in model_order:
        m_tn = m_fp = m_fn = m_tp = 0
        model_thresholds: list[float] = []
        m_rows = [r for r in fold_rows if str(r.get("model_name", "")) == model_name]

        for r in m_rows:
            thr = float(r.get("threshold", 0.5) if r.get("threshold") is not None else 0.5)
            model_thresholds.append(thr)
            yt = np.asarray(r.get("y_true", []), dtype=int)
            yp = np.asarray(r.get("y_pred", []), dtype=float)
            if yt.size == 0 or yp.size == 0:
                continue
            yh = (yp >= thr).astype(int)
            m_tn += int(np.sum((yt == 0) & (yh == 0)))
            m_fp += int(np.sum((yt == 0) & (yh == 1)))
            m_fn += int(np.sum((yt == 1) & (yh == 0)))
            m_tp += int(np.sum((yt == 1) & (yh == 1)))

        confusion_by_model[model_name] = [[m_tn, m_fp], [m_fn, m_tp]]
        if model_thresholds:
            threshold_summary.append(
                {
                    "model_name": model_name,
                    "mean_threshold": float(np.mean(model_thresholds)),
                    "std_threshold": float(np.std(model_thresholds)),
                    "n_folds": int(len(model_thresholds)),
                }
            )

    # Identify the best-performing model for threshold highlighting.
    best_threshold: dict[str, Any] | None = None
    if threshold_summary:
        metric_by_model: dict[str, float] = {}
        for model_name in model_order:
            vals = [
                float(r.get("scores", {}).get(primary_metric, np.nan))
                for r in fold_rows
                if str(r.get("model_name", "")) == model_name
            ]
            vals = [v for v in vals if not np.isnan(v)]
            if vals:
                metric_by_model[model_name] = float(np.mean(vals))

        if metric_by_model:
            lower_is_better = primary_metric in {"rmse", "mad"}
            best_model = (
                min(metric_by_model, key=metric_by_model.get)
                if lower_is_better
                else max(metric_by_model, key=metric_by_model.get)
            )
            best_threshold = next(
                (r for r in threshold_summary if str(r.get("model_name", "")) == best_model),
                None,
            )

    return confusion_by_model, threshold_summary, best_threshold


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def build_task_plot_payload(task_result: dict[str, Any]) -> dict[str, Any]:
    """Build the plotting payload for one task result.

    Handles edge cases (empty folds, all-NaN scores) and produces a
    JSON-serialisable dict consumed by the Plotly JavaScript renderer.
    """
    fold_rows = task_result.get("fold_results", [])
    primary_metric = task_result.get("primary_metric", "metric")
    task = task_result.get("task", "unknown")

    # Edge case: no fold results at all.
    if not fold_rows:
        return _safe_empty_payload(task_result, None, "No fold results available")

    # Build model-fold long-format rows (skip NaN primary scores).
    model_fold_points = []
    for r in fold_rows:
        primary_score = r.get("scores", {}).get(primary_metric, float("nan"))
        if not np.isnan(primary_score):
            model_fold_points.append(
                {
                    "model_name": str(r.get("model_name", "unknown")),
                    "outer_fold": int(r.get("outer_fold", 0)),
                    "primary_score": float(primary_score),
                }
            )

    # Derive display order from first appearance in model_fold_points.
    model_order: list[str] = []
    for row in model_fold_points:
        if row["model_name"] not in model_order:
            model_order.append(row["model_name"])

    # Edge case: all primary scores were NaN.
    if not model_order:
        return _safe_empty_payload(
            task_result, model_fold_points, "All fold scores are NaN"
        )

    # Collect all metric names from fold scores.
    metric_names_all = sorted(
        {str(k) for r in fold_rows for k in r.get("scores", {}).keys() if k is not None}
    )

    # Build metric heatmaps.
    metric_heatmap = _build_metric_heatmap(fold_rows, model_order, metric_names_all)
    metric_heatmap_original = _EMPTY_HEATMAP.copy()
    metric_heatmap_logit = _EMPTY_HEATMAP.copy()

    if task == "regression":
        original_metric_names = sorted(
            [m for m in metric_names_all if m.startswith("original_")]
        )
        logit_metric_names = sorted(
            [m for m in metric_names_all if m.startswith("logit_")]
        )
        # Backward-compatible fallback when prefixed metrics are absent.
        if not original_metric_names:
            original_metric_names = sorted(
                [m for m in metric_names_all if m in {"rmse", "mad", "r2_rss", "ccc"}]
            )
        # Backfill logit metrics for older artifacts that only have original-scale scores.
        if not logit_metric_names:
            for r in fold_rows:
                y_true_row = r.get("y_true", [])
                y_pred_row = r.get("y_pred", [])
                if y_true_row and y_pred_row:
                    y_true_log = _safe_logit(y_true_row)
                    y_pred_log = _safe_logit(y_pred_row)
                    reg_log = regression_metrics(
                        np.asarray(y_true_log, dtype=float),
                        np.asarray(y_pred_log, dtype=float),
                    )
                    s = r.get("scores", {})
                    for k, v in reg_log.items():
                        s[f"logit_{k}"] = float(v)
                    r["scores"] = s
            metric_names_all = sorted(
                {str(k) for r in fold_rows for k in r.get("scores", {}).keys() if k is not None}
            )
            logit_metric_names = sorted(
                [m for m in metric_names_all if m.startswith("logit_")]
            )

        metric_heatmap_original = _build_metric_heatmap(
            fold_rows, model_order, original_metric_names
        )
        metric_heatmap_logit = _build_metric_heatmap(
            fold_rows, model_order, logit_metric_names
        )
        # Suppress the mixed heatmap for regression (use scale-specific ones).
        metric_heatmap = _EMPTY_HEATMAP.copy()

    # Fold-connecting line segments.
    fold_line_rows: dict[int, dict[str, float]] = {}
    for row in model_fold_points:
        fold_line_rows.setdefault(int(row["outer_fold"]), {})[row["model_name"]] = float(
            row["primary_score"]
        )
    fold_lines = [
        {
            "outer_fold": fold_id,
            "x": [m for m in model_order if m in model_to_score],
            "y": [model_to_score[m] for m in model_order if m in model_to_score],
        }
        for fold_id, model_to_score in sorted(fold_line_rows.items())
    ]

    # Prediction and threshold payloads.
    (
        y_true,
        y_pred,
        thresholds,
        prediction_by_model,
        prediction_by_model_scale,
    ) = _build_predictions_payload(fold_rows, model_order, task)

    confusion_by_model: dict[str, list[list[int]]] = {}
    threshold_summary: list[dict[str, Any]] = []
    best_threshold: dict[str, Any] | None = None
    roc_by_model: dict[str, dict[str, Any]] = {}
    if task == "classification":
        confusion_by_model, threshold_summary, best_threshold = (
            _build_classification_threshold_payload(
                fold_rows, model_order, primary_metric
            )
        )
        roc_by_model = _build_roc_payload(fold_rows, model_order)
        pr_by_model = _build_pr_payload(fold_rows, model_order)

    return {
        "task": task,
        "status": task_result.get("status", "unknown"),
        "primary_metric": primary_metric,
        "n_samples": task_result.get("n_samples", None),
        "warnings": task_result.get("warnings", []),
        "model_order": model_order,
        "metric_heatmap": metric_heatmap,
        "metric_heatmap_original": metric_heatmap_original,
        "metric_heatmap_logit": metric_heatmap_logit,
        "model_fold_points": model_fold_points,
        "fold_lines": fold_lines,
        "y_true": y_true,
        "y_pred": y_pred,
        "prediction_by_model": prediction_by_model,
        "prediction_by_model_scale": prediction_by_model_scale,
        "thresholds": thresholds,
        "threshold_by_model": threshold_summary,
        "best_threshold": best_threshold,
        "confusion_matrix": None,
        "confusion_by_model": confusion_by_model,
        "roc_by_model": roc_by_model,
        "pr_by_model": pr_by_model,
        "response_distribution": _ensure_thresholds(task_result),
        "important_params": important_params_table_rows(task_result),
    }
