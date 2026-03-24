from __future__ import annotations

"""Model selection, hyperparameter tuning, and fold evaluation.

This module is the public API shim for the ``splicing_ml.models`` subpackage.
All implementation logic lives in the subpackage; this file:
  - Keeps small utilities that have no natural sub-home
    (split helpers, build_baseline, metric_key, scorer_name).
  - Re-exports the full historical public API so callers and tests that
    import from ``splicing_ml.modeling`` do not need to change.
"""

import math
from typing import Any

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier, DummyRegressor

from .config import RNG_SEED

# ---------------------------------------------------------------------------
# Re-exports from the models subpackage
# ---------------------------------------------------------------------------
from .models.beta import BetaRegressor, _inverse_logit, _logit_transform
from .models.evaluator import evaluate_outer_fold, tune_threshold_balanced_accuracy
from .models.grids import (
    _axis_count_from_budget,
    _lhs_unit,
    _map_with_scale,
    _pick_evenly_spaced,
    build_param_candidates,
    choose_param_grid,
    choose_param_lhs_candidates,
)
from .models.search import fit_best_estimator
from .models.xgb_utils import (
    _set_xgb_cpu_predictor_for_inference,
    _unwrap_model_step,
    xgb_gpu_available,
)

# metric_key and scorer_name are also re-exported for backward compat.
from .models.search import metric_key, scorer_name

__all__ = [
    # Split helpers (defined here)
    "balanced_group_split_indices",
    "bounded_group_splits",
    "inner_split_indices",
    # Baseline builder (defined here)
    "build_baseline",
    # Metric helpers
    "metric_key",
    "scorer_name",
    # Re-exports from models subpackage
    "BetaRegressor",
    "_inverse_logit",
    "_logit_transform",
    "evaluate_outer_fold",
    "tune_threshold_balanced_accuracy",
    "build_param_candidates",
    "choose_param_grid",
    "choose_param_lhs_candidates",
    "fit_best_estimator",
    "_set_xgb_cpu_predictor_for_inference",
    "_unwrap_model_step",
    "xgb_gpu_available",
    # Grid helpers
    "_axis_count_from_budget",
    "_lhs_unit",
    "_map_with_scale",
    "_pick_evenly_spaced",
]


# ---------------------------------------------------------------------------
# Split helpers — kept here because they are tested via splicing_ml.modeling
# and are also used directly by pipeline.py
# ---------------------------------------------------------------------------


def bounded_group_splits(
    requested_splits: int, n_groups: int, max_splits: int = 10
) -> int:
    """Clamp grouped CV fold count to [2, max_splits] and available groups."""
    return max(2, min(requested_splits, n_groups, max_splits))


def balanced_group_split_indices(
    groups: pd.Series,
    n_splits: int,
    seed: int = RNG_SEED,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Create grouped CV splits balanced by number of observations.

    Groups are assigned to folds greedily by descending group size
    (bin-packing) so each fold has approximately balanced observation counts.
    Group exclusivity is enforced: an entire group is always in either train
    or test, never both.
    """
    g = groups.reset_index(drop=True)
    unique_groups = g.dropna().unique().tolist()
    if len(unique_groups) < 2:
        raise ValueError("Need at least 2 unique groups for grouped CV")

    n_splits = bounded_group_splits(n_splits, len(unique_groups), max_splits=10)
    group_counts = g.value_counts().to_dict()

    rng = np.random.default_rng(seed)
    shuffled = list(unique_groups)
    rng.shuffle(shuffled)
    # Sort descending by size for the bin-packing pass.
    shuffled.sort(key=lambda k: group_counts.get(k, 0), reverse=True)

    fold_groups: list[list[Any]] = [[] for _ in range(n_splits)]
    fold_sizes = [0 for _ in range(n_splits)]

    # Greedy bin-packing: assign each group to the fold with the fewest rows.
    for grp in shuffled:
        target_fold = int(np.argmin(fold_sizes))
        fold_groups[target_fold].append(grp)
        fold_sizes[target_fold] += int(group_counts.get(grp, 0))

    all_idx = np.arange(g.shape[0])
    splits: list[tuple[np.ndarray, np.ndarray]] = []
    for test_groups in fold_groups:
        if not test_groups:
            continue
        test_mask = g.isin(test_groups).to_numpy()
        test_idx = all_idx[test_mask]
        train_idx = all_idx[~test_mask]
        if train_idx.size == 0 or test_idx.size == 0:
            continue
        splits.append((train_idx, test_idx))

    if len(splits) < 2:
        raise ValueError("Could not build at least 2 non-empty grouped folds")
    return splits


def inner_split_indices(
    groups: pd.Series, n_splits: int
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Build observation-balanced grouped splits for inner CV.

    Thin wrapper around ``balanced_group_split_indices``.
    """
    return balanced_group_split_indices(groups, n_splits=n_splits)


# ---------------------------------------------------------------------------
# Baseline builder
# ---------------------------------------------------------------------------


def build_baseline(task: str):
    """Return a task-appropriate baseline estimator for delta reporting."""
    return (
        DummyClassifier(strategy="prior")
        if task == "classification"
        else DummyRegressor(strategy="mean")
    )


# ---------------------------------------------------------------------------
# Legacy search-strategy helper — kept for backward compat
# ---------------------------------------------------------------------------


def _use_random_search(search_strategy: str, model_name: str) -> bool:
    """Return True when RandomizedSearchCV / LHS should be used."""
    if search_strategy == "random":
        return True
    if search_strategy == "grid":
        return False
    return model_name in {"rf", "xgb"}


# ---------------------------------------------------------------------------
# Legacy helpers still referenced in some contexts
# ---------------------------------------------------------------------------


def _model_outputs_psi_scale(estimator: Any) -> bool:
    """Return True when the estimator predicts on the PSI [0, 1] scale."""
    from .models.evaluator import _model_outputs_psi_scale as _impl

    return _impl(estimator)
