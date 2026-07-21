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


MAX_GLOBAL_POINTS = 8000
MAX_PER_MODEL_POINTS = 2000
MAX_PSI_SAMPLE_POINTS = 30000
MAX_CURVE_POINTS = 3000
MAX_FEATURE_DIST_POINTS = 3000


def _downsample_pair(
    x_vals: list[float], y_vals: list[float], max_points: int
) -> tuple[list[float], list[float]]:
    """Downsample paired vectors using evenly spaced deterministic indices."""
    n = min(len(x_vals), len(y_vals))
    if n <= max_points:
        return x_vals[:n], y_vals[:n]
    idx = np.linspace(0, n - 1, num=max_points, dtype=int)
    return [x_vals[i] for i in idx], [y_vals[i] for i in idx]


def _downsample_list(vals: list[float], max_points: int) -> list[float]:
    """Downsample one vector using evenly spaced deterministic indices."""
    if len(vals) <= max_points:
        return vals
    idx = np.linspace(0, len(vals) - 1, num=max_points, dtype=int)
    return [vals[i] for i in idx]


def _downsample_curve(
    x_vals: list[float], y_vals: list[float], max_points: int
) -> tuple[list[float], list[float]]:
    """Downsample curve points while preserving paired x/y ordering."""
    return _downsample_pair(x_vals, y_vals, max_points)


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

_LOWER_IS_BETTER_METRICS = {"rmse", "mad"}


def _is_lower_is_better(metric: str) -> bool:
    """True if a smaller value is better for *metric* (e.g. RMSE, MAD).

    Handles the ``original_``/``logit_`` scale prefixes used for regression
    metrics; everything else (r2_rss, ccc, and all classification metrics)
    is higher-is-better.
    """
    base = metric
    for prefix in ("original_", "logit_"):
        if base.startswith(prefix):
            base = base[len(prefix) :]
            break
    return base in _LOWER_IS_BETTER_METRICS


def _direction_sort_key(metric: str) -> tuple[int, str]:
    """Sort key grouping lower-is-better metrics before higher-is-better ones."""
    return (0 if _is_lower_is_better(metric) else 1, metric)


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
        "directions": [_is_lower_is_better(m) for m in metric_names],
    }


_EMPTY_HEATMAP: dict[str, Any] = {
    "models": [],
    "metrics": [],
    "z": [],
    "sd": [],
    "text": [],
    "directions": [],
}


def _build_metric_fold_data(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
    metric_names: list[str],
) -> dict[str, dict[str, Any]]:
    """Build per-metric fold-level box-plot points and fold-connecting lines.

    Mirrors the primary-metric-only ``model_fold_points``/``fold_lines``
    computation, but repeated for every metric so the HTML report can render
    one box plot per metric instead of only the primary one.
    """
    out: dict[str, dict[str, Any]] = {}
    for metric in metric_names:
        points: list[dict[str, Any]] = []
        for r in fold_rows:
            score = float(r.get("scores", {}).get(metric, float("nan")))
            if not np.isnan(score):
                points.append(
                    {
                        "model_name": str(r.get("model_name", "unknown")),
                        "outer_fold": int(r.get("outer_fold", 0)),
                        "score": score,
                    }
                )
        fold_to_scores: dict[int, dict[str, float]] = {}
        for p in points:
            fold_to_scores.setdefault(p["outer_fold"], {})[p["model_name"]] = p[
                "score"
            ]
        lines = [
            {
                "outer_fold": fold_id,
                "x": [m for m in model_order if m in model_to_score],
                "y": [model_to_score[m] for m in model_order if m in model_to_score],
            }
            for fold_id, model_to_score in sorted(fold_to_scores.items())
        ]
        out[metric] = {"points": points, "lines": lines}
    return out


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
        "metric_fold_data": {},
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
        "pr_prevalence": None,
        "response_distribution": _ensure_thresholds(task_result),
        "transformed_feature_distributions": {},
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
    psi_sample = dist.get("psi_sample")
    if isinstance(psi_sample, list):
        dist["psi_sample"] = _downsample_list(psi_sample, MAX_PSI_SAMPLE_POINTS)
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
            if len(row_y_true_logit) == 0 and row_y_true:
                row_y_true_logit = _safe_logit(row_y_true)
            if len(row_y_pred_logit) == 0 and row_y_pred:
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
            prediction_by_model_scale[row_model]["original"]["y_true"].extend(
                row_y_true
            )
            prediction_by_model_scale[row_model]["original"]["y_pred"].extend(
                row_y_pred
            )
            prediction_by_model_scale[row_model]["logit"]["y_true"].extend(
                row_y_true_logit
            )
            prediction_by_model_scale[row_model]["logit"]["y_pred"].extend(
                row_y_pred_logit
            )

    # Down-sample to keep HTML payload lightweight for large runs.
    y_true, y_pred = _downsample_pair(y_true, y_pred, MAX_GLOBAL_POINTS)

    for model_name in model_order:
        model_pred = prediction_by_model.get(model_name, {})
        ds_true, ds_pred = _downsample_pair(
            model_pred.get("y_true", []),
            model_pred.get("y_pred", []),
            MAX_PER_MODEL_POINTS,
        )
        model_pred["y_true"] = ds_true
        model_pred["y_pred"] = ds_pred
        prediction_by_model[model_name] = model_pred

        model_scales = prediction_by_model_scale.get(model_name, {})
        for scale_name in ("original", "logit"):
            scale_vals = model_scales.get(scale_name, {"y_true": [], "y_pred": []})
            s_true, s_pred = _downsample_pair(
                scale_vals.get("y_true", []),
                scale_vals.get("y_pred", []),
                MAX_PER_MODEL_POINTS,
            )
            scale_vals["y_true"] = s_true
            scale_vals["y_pred"] = s_pred
            model_scales[scale_name] = scale_vals
        prediction_by_model_scale[model_name] = model_scales

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
        fpr_list, tpr_list = _downsample_curve(
            fpr.tolist(), tpr.tolist(), MAX_CURVE_POINTS
        )
        auc_val = float(auc(fpr, tpr))
        roc_by_model[model_name] = {
            "fpr": fpr_list,
            "tpr": tpr_list,
            "auc": auc_val,
        }
    return roc_by_model


def _build_pr_payload(
    fold_rows: list[dict[str, Any]],
    model_order: list[str],
) -> tuple[dict[str, dict[str, Any]], float | None]:
    """Build per-model PR curve data by pooling predictions across folds.

    Returns a tuple of:
    - dict keyed by model name, each value containing ``precision``,
      ``recall`` (lists), and ``avg_precision`` (AUPR float).
      Includes a "Baseline (prior)" entry when ``baseline_y_pred`` is available.
    - prevalence (float) of the positive class across all folds, or None.

    Models with fewer than two classes in the pooled labels are skipped.
    """
    from sklearn.metrics import auc, precision_recall_curve

    pr_by_model: dict[str, dict[str, Any]] = {}
    all_y_true_pooled: list[float] = []

    for model_name in model_order:
        m_rows = [r for r in fold_rows if str(r.get("model_name", "")) == model_name]
        all_y_true: list[float] = []
        all_y_pred: list[float] = []
        for r in m_rows:
            all_y_true.extend(r.get("y_true", []))
            all_y_pred.extend(r.get("y_pred", []))
        if not all_y_true:
            continue
        if not all_y_true_pooled:
            all_y_true_pooled = all_y_true
        yt = np.asarray(all_y_true, dtype=int)
        yp = np.asarray(all_y_pred, dtype=float)
        if len(np.unique(yt)) < 2:
            continue
        precision, recall, _ = precision_recall_curve(yt, yp)
        recall_list, precision_list = _downsample_curve(
            recall.tolist(), precision.tolist(), MAX_CURVE_POINTS
        )
        aupr = float(auc(recall, precision))
        pr_by_model[model_name] = {
            "precision": precision_list,
            "recall": recall_list,
            "avg_precision": aupr,
        }

    prevalence: float | None = (
        float(np.mean(np.asarray(all_y_true_pooled, dtype=int)))
        if all_y_true_pooled
        else None
    )

    # Baseline model (DummyClassifier prior) PR curve.
    baseline_yt: list[float] = []
    baseline_yp: list[float] = []
    for r in fold_rows:
        bpred = r.get("baseline_y_pred", [])
        yt_row = r.get("y_true", [])
        if bpred and yt_row:
            baseline_yt.extend(yt_row)
            baseline_yp.extend(bpred)
    if baseline_yt and len(np.unique(np.asarray(baseline_yt, dtype=int))) >= 2:
        byt = np.asarray(baseline_yt, dtype=int)
        byp = np.asarray(baseline_yp, dtype=float)
        precision, recall, _ = precision_recall_curve(byt, byp)
        recall_list, precision_list = _downsample_curve(
            recall.tolist(), precision.tolist(), MAX_CURVE_POINTS
        )
        aupr = float(auc(recall, precision))
        pr_by_model["Baseline (prior)"] = {
            "precision": precision_list,
            "recall": recall_list,
            "avg_precision": aupr,
        }

    return pr_by_model, prevalence


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
            thr = float(
                r.get("threshold", 0.5) if r.get("threshold") is not None else 0.5
            )
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
                (
                    r
                    for r in threshold_summary
                    if str(r.get("model_name", "")) == best_model
                ),
                None,
            )

    return confusion_by_model, threshold_summary, best_threshold


def _build_transformed_feature_distributions(
    fold_rows: list[dict[str, Any]],
) -> dict[str, list[float]]:
    """Pool transformed-feature samples across every fold/model into one
    per-feature list, for a post-preprocessing distribution histogram.

    The preprocessor is refit per fold/model, but on the same fold's
    training data with the same transform steps -- pooling across models is
    harmless (near-duplicate samples of the same underlying population) and
    gives a bigger, smoother sample than picking just one model.
    """
    pooled: dict[str, list[float]] = {}
    for r in fold_rows:
        tf = r.get("transformed_features")
        if not tf:
            continue
        names = tf.get("feature_names") or []
        values = tf.get("values") or []
        if not names or not values:
            continue
        arr = np.asarray(values, dtype=float)
        for i, name in enumerate(names):
            pooled.setdefault(name, []).extend(arr[:, i].tolist())
    return {
        name: _downsample_list(vals, MAX_FEATURE_DIST_POINTS)
        for name, vals in pooled.items()
    }


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

    # Inject synthetic baseline fold rows so the baseline appears in the
    # metric heatmap and scatter plots. One row per outer fold, derived from
    # the first available model row for that fold.
    _baseline_label = "Baseline (mean)" if task == "regression" else "Baseline (prior)"
    _seen_folds: set[int] = set()
    for r in list(fold_rows):
        fid = int(r.get("outer_fold", -1))
        if fid in _seen_folds:
            continue
        bs = r.get("baseline_scores", {})
        if not bs:
            continue
        _seen_folds.add(fid)
        fold_rows = list(fold_rows) + [
            {
                "model_name": _baseline_label,
                "scores": bs,
                "outer_fold": fid,
                "y_true": r.get("y_true", []),
                "y_pred": r.get("baseline_y_pred", []),
            }
        ]

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

    # Build metric heatmaps. Metrics are grouped lower-is-better-first so the
    # heatmap/box-plot grid visually separates error metrics from fit-quality
    # ones instead of interleaving directions the row-color scale can't both
    # satisfy at once.
    metric_names_all = sorted(metric_names_all, key=_direction_sort_key)
    metric_heatmap = _build_metric_heatmap(fold_rows, model_order, metric_names_all)
    metric_heatmap_original = _EMPTY_HEATMAP.copy()
    metric_heatmap_logit = _EMPTY_HEATMAP.copy()

    if task == "regression":
        original_metric_names = sorted(
            [m for m in metric_names_all if m.startswith("original_")],
            key=_direction_sort_key,
        )
        logit_metric_names = sorted(
            [m for m in metric_names_all if m.startswith("logit_")],
            key=_direction_sort_key,
        )
        # Backward-compatible fallback when prefixed metrics are absent.
        if not original_metric_names:
            original_metric_names = sorted(
                [m for m in metric_names_all if m in {"rmse", "mad", "r2_rss", "ccc"}],
                key=_direction_sort_key,
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
                {
                    str(k)
                    for r in fold_rows
                    for k in r.get("scores", {}).keys()
                    if k is not None
                }
            )
            logit_metric_names = sorted(
                [m for m in metric_names_all if m.startswith("logit_")],
                key=_direction_sort_key,
            )

        metric_heatmap_original = _build_metric_heatmap(
            fold_rows, model_order, original_metric_names
        )
        metric_heatmap_logit = _build_metric_heatmap(
            fold_rows, model_order, logit_metric_names
        )
        # Suppress the mixed heatmap for regression (use scale-specific ones).
        metric_heatmap = _EMPTY_HEATMAP.copy()
        metric_fold_data = _build_metric_fold_data(
            fold_rows, model_order, original_metric_names + logit_metric_names
        )
    else:
        metric_fold_data = _build_metric_fold_data(
            fold_rows, model_order, metric_names_all
        )

    # Fold-connecting line segments.
    fold_line_rows: dict[int, dict[str, float]] = {}
    for row in model_fold_points:
        fold_line_rows.setdefault(int(row["outer_fold"]), {})[row["model_name"]] = (
            float(row["primary_score"])
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
    pr_by_model: dict[str, dict[str, Any]] = {}
    pr_prevalence: float | None = None
    if task == "classification":
        confusion_by_model, threshold_summary, best_threshold = (
            _build_classification_threshold_payload(
                fold_rows, model_order, primary_metric
            )
        )
        roc_by_model = _build_roc_payload(fold_rows, model_order)
        pr_by_model, pr_prevalence = _build_pr_payload(fold_rows, model_order)

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
        "metric_fold_data": metric_fold_data,
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
        "pr_prevalence": pr_prevalence,
        "response_distribution": _ensure_thresholds(task_result),
        "transformed_feature_distributions": _build_transformed_feature_distributions(
            fold_rows
        ),
        "important_params": important_params_table_rows(task_result),
    }
