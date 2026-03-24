from __future__ import annotations

"""Metric utilities for splicing ML tasks.

This module centralizes scoring logic for both regression and classification
workflows so all folds/configurations are evaluated consistently.
"""

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    mean_squared_error,
    roc_auc_score,
)


__all__ = [
    "r2_from_rss",
    "concordance_correlation_coefficient",
    "regression_metrics",
    "classification_metrics",
    "bootstrap_ci",
]


def r2_from_rss(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute RSS-based R-squared.

    This implementation uses the residual sum of squares (RSS) definition,
    not Pearson correlation squared, matching the requested specification.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    rss = np.sum((y_true - y_pred) ** 2)
    tss = np.sum((y_true - np.mean(y_true)) ** 2)
    if tss == 0:
        return 0.0
    return 1.0 - (rss / tss)


def concordance_correlation_coefficient(
    y_true: np.ndarray, y_pred: np.ndarray
) -> float:
    """Compute Lin's concordance correlation coefficient (CCC).

    CCC measures both correlation and agreement around the identity line.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)
    cov = np.mean((y_true - mean_true) * (y_pred - mean_pred))
    denom = var_true + var_pred + (mean_true - mean_pred) ** 2
    if denom == 0:
        return 0.0
    return float((2.0 * cov) / denom)


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Return all requested regression metrics in a single dictionary."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    mad = float(np.mean(np.abs(y_true - y_pred)))
    r2_rss = float(r2_from_rss(y_true, y_pred))
    ccc = float(concordance_correlation_coefficient(y_true, y_pred))
    return {
        "rmse": rmse,
        "mad": mad,
        "r2_rss": r2_rss,
        "ccc": ccc,
    }


def classification_metrics(
    y_true: np.ndarray, y_prob: np.ndarray, threshold: float
) -> dict[str, float]:
    """Return requested classification metrics from probabilities and threshold.

    Parameters
    ----------
    y_true
        Binary ground truth labels (0/1).
    y_prob
        Positive-class probabilities.
    threshold
        Decision threshold used to derive hard labels.
    """
    y_true = np.asarray(y_true, dtype=int)
    y_prob = np.asarray(y_prob, dtype=float)
    y_hat = (y_prob >= threshold).astype(int)

    metrics: dict[str, float] = {
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_hat)),
        "auprc": float(average_precision_score(y_true, y_prob)),
        "f1": float(f1_score(y_true, y_hat, zero_division=0)),
        "mcc": float(matthews_corrcoef(y_true, y_hat)),
    }

    # AUROC requires both classes to be present in y_true.
    if len(np.unique(y_true)) == 2:
        metrics["auroc"] = float(roc_auc_score(y_true, y_prob))
    else:
        metrics["auroc"] = float("nan")

    return metrics


def bootstrap_ci(
    values: np.ndarray,
    n_bootstrap: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> tuple[float, float]:
    """Compute percentile bootstrap confidence interval for the mean.

    This is used on outer-fold primary metric values to report uncertainty.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    boot_means = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        sample = rng.choice(values, size=values.size, replace=True)
        boot_means[i] = np.mean(sample)
    low = float(np.quantile(boot_means, alpha / 2.0))
    high = float(np.quantile(boot_means, 1.0 - alpha / 2.0))
    return (low, high)
