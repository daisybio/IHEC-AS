from __future__ import annotations

"""Hyperparameter tuning and best-estimator selection.

Public entry point: ``fit_best_estimator`` — dispatches to one of three
strategy helpers based on model type:

- ``_fit_elasticnet_cv``  : ElasticNetCV / LogisticRegressionCV self-select
                            alpha (or C) and l1_ratio via built-in CV.
- ``_fit_beta``           : BetaRegressor has no hyperparameters; uses
                            ``sklearn.model_selection.cross_validate`` to score
                            then refits on the full training fold.
- ``_fit_grid_search``    : All other models use GridSearchCV (or LHS-based
                            GridSearchCV) with optional GPU retry and XGBoost
                            post-search refit.
"""

import json
import math
import os
from typing import Any
import warnings

import numpy as np
import pandas as pd
from joblib import parallel_backend
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import ElasticNetCV
from sklearn.exceptions import FitFailedWarning
from sklearn.model_selection import GridSearchCV, ParameterGrid, cross_validate
from sklearn.pipeline import Pipeline

from ..config import RNG_SEED
from ..utils import progress_iter, safe_json, sanitize_best_params, vlog
from .beta import BetaRegressor, _inverse_logit
from ._contexts import _es_raw_val_context, _es_pp_val_context
from .grids import (
    build_param_candidates,
    _use_lhs_strategy,
    _ALPHA_GRID,
    _C_GRID,
    _L1_RATIO_GRID,
    MLP_ES_PATIENCE,
)
from .xgb_utils import _set_xgb_cpu_predictor_for_inference, xgb_gpu_available
from .lgbm_utils import _set_lgbm_cpu_predictor_for_inference, lgbm_gpu_available

__all__ = ["fit_best_estimator"]


class _ProgressGridSearchCV(GridSearchCV):
    """GridSearchCV variant that can show candidate progress in debug runs.

    When enabled, candidates are evaluated one-by-one so a progress bar can be
    updated after each completed candidate. This is intended for debug
    observability and is disabled by default.
    """

    def __init__(
        self,
        estimator,
        param_grid,
        *,
        scoring=None,
        n_jobs=None,
        refit=True,
        cv=None,
        verbose=0,
        pre_dispatch="2*n_jobs",
        error_score=np.nan,
        return_train_score=False,
        progress_enabled: bool = False,
        progress_desc: str | None = None,
    ) -> None:
        """Initialize a _ProgressGridSearchCV instance."""
        super().__init__(
            estimator=estimator,
            param_grid=param_grid,
            scoring=scoring,
            n_jobs=n_jobs,
            refit=refit,
            cv=cv,
            verbose=verbose,
            pre_dispatch=pre_dispatch,
            error_score=error_score,
            return_train_score=return_train_score,
        )
        # Keep sklearn BaseEstimator param introspection happy by storing
        # constructor params under identical attribute names.
        self.progress_enabled = bool(progress_enabled)
        self.progress_desc = progress_desc

    def _run_search(self, evaluate_candidates) -> None:  # type: ignore[override]
        """Internal helper for run search."""
        candidates = list(ParameterGrid(self.param_grid))
        if not self.progress_enabled:
            evaluate_candidates(candidates)
            return

        for candidate in progress_iter(
            candidates,
            total=len(candidates),
            desc=self.progress_desc,
            enabled=True,
        ):
            evaluate_candidates([candidate])


def scorer_name(task: str) -> str:
    """Return the sklearn scorer name aligned with the primary optimization metric."""
    # This is consumed by sklearn search APIs (GridSearchCV/cross_validate/OptunaSearchCV),
    # so we must return sklearn scorer identifiers, including sign conventions.
    return "roc_auc" if task == "classification" else "neg_root_mean_squared_error"


def metric_key(task: str) -> str:
    """Return the primary metric name for a task."""
    # This key is used for internal result dictionaries/logging, not sklearn scorer names.
    return "auroc" if task == "classification" else "rmse"


def _adaptive_c_grid(n_samples: int) -> tuple[list[float], list[float]]:
    """Return (C_grid, l1_ratio_grid) scaled to the training set size.

    Fewer candidates for larger n_samples to maintain practical runtimes.
    The C grid is always log-spaced over [1e-6, 1e2]; only its density changes.
    l1_ratio grid is thinned for very large datasets.
    """
    if n_samples < 50_000:
        n_cs, l1_ratios = 100, [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]
    elif n_samples < 200_000:
        n_cs, l1_ratios = 50, [0.1, 0.5, 0.9, 0.95, 1.0]
    elif n_samples < 600_000:
        n_cs, l1_ratios = 25, [0.1, 0.5, 0.9, 1.0]
    else:
        n_cs, l1_ratios = 25, [0.1, 0.9, 1.0]
    c_grid = (1.0 / np.logspace(-6, 2, num=n_cs)).tolist()
    return c_grid, l1_ratios


# ---------------------------------------------------------------------------
# Private strategy helpers
# ---------------------------------------------------------------------------


def _fit_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    verbose: bool,
    cuml_use_gpu: bool | None = None,
    debug_grid_progress: bool = False,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit ElasticNetCV — alpha and l1_ratio selection via built-in CV.

    The same grouped inner-fold splits used by all other models are passed
    directly to ElasticNetCV so fold assignments are consistent across models.
    Always sweeps _ALPHA_GRID alphas and _L1_RATIO_GRID l1_ratios.

    The cuML GPU path (_fit_cuml_elasticnet_cv) is intentionally bypassed:
    ElasticNetCV uses warm-start coordinate descent, so sweeping 100 alphas costs
    nearly the same as a small slice (one solver pass per l1_ratio with warm
    starts).  Benchmarks on 477k samples showed CPU (32 cores) at ~20 min vs
    GPU at ~30 min for 105 cold cuML fits, with GPU also scoring ~3% lower due
    to float32 precision and max_iter=5000 non-convergence.
    At 16 cores the wall times converge (~40 min CPU vs ~30 min GPU), but the
    accuracy gap remains — so CPU is preferred unconditionally until a
    regularization-path cuML variant becomes available.
    """
    # ElasticNetCV uses warm-start coordinate descent — sweeping the
    # full alpha/l1_ratio grid costs nearly the same as a tiny slice.
    c_grid, l1_grid = _adaptive_c_grid(x_train.shape[0])
    # ElasticNetCV uses alpha, which is the reciprocal of C
    alpha_grid = [float(1.0 / c) for c in c_grid]
    l1_grid = [float(r) for r in l1_grid]

    # ElasticNetCV parallelises over n_l1_ratios tasks (each task sweeps all
    # alphas for all cv folds via warm-start coordinate descent).
    # Empirical reference: ~10 min/batch at 477k samples and 100 alphas.
    _n_batches = math.ceil(len(l1_grid) / max(1, max_cores))
    _est_min = (
        _n_batches * 10.0 * (x_train.shape[0] / 477_000) * (len(alpha_grid) / 100)
    )
    vlog(
        verbose,
        f"Fitting elasticnet (ElasticNetCV), task={task}, "
        f"n_samples={x_train.shape[0]}, n_features={x_train.shape[1]}, "
        f"n_alphas={len(alpha_grid)}, n_l1_ratios={len(l1_grid)}, "
        f"batches={_n_batches} (cores={max_cores}), ~{_est_min:.0f} min est.",
    )
    estimator = Pipeline(
        steps=[
            ("prep", preprocessor),
            (
                "model",
                ElasticNetCV(
                    random_state=RNG_SEED,
                    max_iter=20000,
                    cv=inner_cv,
                    alphas=alpha_grid,
                    l1_ratio=l1_grid,
                    n_jobs=max_cores,
                    verbose=1 if debug_grid_progress else 0,
                ),
            ),
        ]
    )
    # Non-convergence at extreme alpha values is expected during the warm-start
    # path sweep and does not affect the selected optimum — suppress the noise.
    from sklearn.exceptions import ConvergenceWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        estimator.fit(x_train, y_train)
    model = estimator.named_steps["model"]

    # Extract cross-validated RMSE from the MSE path produced by ElasticNetCV.
    mse_path = np.asarray(getattr(model, "mse_path_", np.array([])))
    best_score = float("nan")
    if mse_path.size:
        best_mse = float(np.min(np.mean(mse_path, axis=1)))
        best_score = -float(np.sqrt(best_mse))

    alpha = float(getattr(model, "alpha_", np.nan))
    l1_ratio = float(getattr(model, "l1_ratio_", np.nan))
    n_alphas = int(np.asarray(getattr(model, "alphas_", np.array([]))).size)
    vlog(
        verbose,
        f"Search complete model=elasticnet, task={task}, strategy=elasticnetcv, "
        f"candidates={n_alphas if n_alphas > 0 else len(_ALPHA_GRID)}, best_score={best_score:.6f}",
    )

    # Replace ElasticNetCV with a plain ElasticNet using the selected alpha and
    # l1_ratio. ElasticNetCV.fit() always re-runs the full warm-start path sweep,
    # so leaving it in place causes the evaluator's refit call to duplicate the
    # search cost. ElasticNet.fit() is a single coordinate descent pass.
    # Use cuML ElasticNet for the final fit when GPU is available — a single
    # fixed-param fit avoids the float32/convergence issues that made the full
    # CV sweep prefer CPU.
    from .cuml_utils import cuml_gpu_available

    if cuml_use_gpu is not False and cuml_gpu_available():
        from cuml.linear_model import ElasticNet as _ElasticNet

        fixed_model = _ElasticNet(
            alpha=alpha, l1_ratio=l1_ratio, max_iter=20000, verbose=False
        )
    else:
        from sklearn.linear_model import ElasticNet as _ElasticNet

        fixed_model = _ElasticNet(
            alpha=alpha, l1_ratio=l1_ratio, max_iter=20000, random_state=RNG_SEED
        )
    estimator = Pipeline(
        steps=[("prep", estimator.named_steps["prep"]), ("model", fixed_model)]
    )

    return estimator, {
        "best_params": sanitize_best_params(
            {
                "model": fixed_model,
                "model__alpha": alpha,
                "model__l1_ratio": l1_ratio,
                "model__cv_folds": int(len(inner_cv)),
                "model__alpha_count": len(alpha_grid),
            }
        ),
        "best_score": best_score,
    }


def _fit_cuml_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    alpha_grid: list[float],
    l1_grid: list[float],
    verbose: bool,
    debug_grid_progress: bool = False,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit cuML ElasticNet via GridSearchCV over a coarser alpha/l1_ratio grid.

    cuML ElasticNet has the same alpha and l1_ratio parameters as sklearn but
    no CV-path variant.  GridSearchCV drives the search (n_jobs=1 for GPU).

    DEPRECATED: not called by _fit_elasticnet_cv.  sklearn ElasticNetCV is faster
    and more accurate (warm-start coordinate descent vs cold cuML float32 fits).
    Retained for reference until a regularization-path cuML variant exists.
    """
    from .cuml_utils import get_cuml_elasticnet

    pipe = Pipeline(
        steps=[
            ("prep", preprocessor),
            ("model", get_cuml_elasticnet()),
        ]
    )
    candidates = [
        {"model__alpha": [float(a)], "model__l1_ratio": [float(r)]}
        for a in alpha_grid
        for r in l1_grid
    ]
    search = _ProgressGridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=inner_cv,
        n_jobs=1,
        refit=True,
        error_score="raise",
        progress_enabled=debug_grid_progress,
        progress_desc=f"elasticnet [{task}]",
    )
    try:
        search.fit(x_train, y_train)
    except Exception as exc:
        if _looks_like_gpu_failure(exc):
            vlog(
                verbose,
                f"cuML ElasticNet GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            # Fall back to the CPU ElasticNetCV path.
            from sklearn.linear_model import ElasticNetCV as _ElasticNetCV

            cpu_est = Pipeline(
                steps=[
                    ("prep", preprocessor),
                    (
                        "model",
                        _ElasticNetCV(
                            random_state=RNG_SEED,
                            max_iter=20000,
                            cv=inner_cv,
                            alphas=alpha_grid,
                            l1_ratio=l1_grid,
                            n_jobs=1,
                        ),
                    ),
                ]
            )
            cpu_est.fit(x_train, y_train)
            cpu_model = cpu_est.named_steps["model"]
            mse_path = np.asarray(getattr(cpu_model, "mse_path_", np.array([])))
            best_score = float("nan")
            if mse_path.size:
                best_mse = float(np.min(np.mean(mse_path, axis=1)))
                best_score = -float(np.sqrt(best_mse))
            alpha = float(getattr(cpu_model, "alpha_", np.nan))
            l1_ratio = float(getattr(cpu_model, "l1_ratio_", np.nan))
            vlog(
                verbose,
                f"Search complete model=elasticnet, task={task}, strategy=elasticnetcv_cpu_fallback, "
                f"candidates={len(alpha_grid) * len(l1_grid)}, best_score={best_score:.6f}",
            )
            return cpu_est, {
                "best_params": sanitize_best_params(
                    {
                        "model": cpu_model,
                        "model__alpha": alpha,
                        "model__l1_ratio": l1_ratio,
                        "model__cv_folds": int(len(inner_cv)),
                        "model__alpha_count": len(alpha_grid),
                    }
                ),
                "best_score": best_score,
            }
        raise

    best_params = search.best_params_
    total_candidates = len(alpha_grid) * len(l1_grid)
    vlog(
        verbose,
        f"Search complete model=elasticnet, task={task}, strategy=cuml_elasticnet_grid, "
        f"candidates={total_candidates}, best_score={float(search.best_score_):.6f}",
    )
    return search.best_estimator_, {
        "best_params": sanitize_best_params(best_params),
        "best_score": float(search.best_score_),
    }


def _fit_logistic_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    verbose: bool,
    cuml_use_gpu: bool | None = None,
    debug_grid_progress: bool = False,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit LogisticRegressionCV — C and l1_ratio selection via built-in CV.

    Mirrors the regression ElasticNetCV path: the same grouped inner-fold splits
    are passed directly so fold assignments are consistent across models, and
    n_jobs is set to max_cores.  Sweeps _C_GRID C values and _L1_RATIO_GRID
    l1_ratios (C = 1 / alpha, reciprocal of _ALPHA_GRID).

    The cuML GPU path (_fit_cuml_logistic_elasticnet_cv) is intentionally bypassed:
    LogisticRegressionCV parallelises over folds × l1_ratios (e.g. 5 × 7 = 35
    tasks) and sweeps all Cs per task via warm-start.  Benchmarks on 477k samples:
    CPU (32 cores) ~20 min / score 0.680 vs GPU 105 cold cuML fits ~30 min /
    score 0.650.  The ~3% accuracy gap is driven by cuML's float32 default and
    max_iter=5000 non-convergence.
    At 16 cores CPU wall time rises to ~40 min (2 batches of 35/16 tasks) vs
    GPU ~30 min, so the crossover is around 20-24 cores — but the accuracy gap
    persists regardless of core count.  CPU is preferred unconditionally until a
    regularization-path cuML variant (analogous to LogisticRegressionCV) exists.
    """
    # LogisticRegressionCV uses warm-start path — sweeping the full
    # C/l1_ratio grid costs nearly the same as a tiny slice.
    c_grid, l1_grid = _adaptive_c_grid(x_train.shape[0])
    c_grid = [float(c) for c in c_grid]
    l1_grid = [float(r) for r in l1_grid]

    from sklearn.linear_model import LogisticRegressionCV

    # Estimate wall time: LogisticRegressionCV parallelises over
    # n_l1_ratios × inner_folds tasks.  Six-point empirical calibration on
    # shared-gpu partition (CPU cores on GPU nodes are slower than CPU nodes):
    #   12k samples,   50 Cs, 2 batches -> ~1.75 min  => ~0.9 min/batch @50Cs
    #   353k samples,  50 Cs, 2 batches -> ~19 min    => ~9.5 min/batch @50Cs
    #   1567k samples, 50 Cs, 1 batch  -> ~63 min     => ~63 min/batch @50Cs
    #   1572k samples, 50 Cs, 1 batch  -> ~79 min     => ~79 min/batch @50Cs
    #   1567k samples, 50 Cs, 1 batch  -> ~62 min     => ~62 min/batch @50Cs
    #   1572k samples, 50 Cs, 1 batch  -> ~83 min     => ~83 min/batch @50Cs
    # Large-run spread (~62-83 min near 1.57M samples) indicates cluster-load
    # variance, so keep a mildly conservative central estimate.
    # Model: (10 + n_samples/11_000) * (n_Cs/100) min/batch.
    # NOTE: For very large datasets (>600k), SAGA solver overhead causes real times to
    # be 1.5-2.3× the estimate. Multiply estimate by 1.8 for best-effort accuracy.
    _n_cv_tasks = len(l1_grid) * len(inner_cv)
    _n_batches = math.ceil(_n_cv_tasks / max(1, max_cores))
    _est_min = _n_batches * (10.0 + x_train.shape[0] / 11_000) * (len(c_grid) / 100)
    vlog(
        verbose,
        f"Fitting elasticnet (LogisticRegressionCV), task={task}, "
        f"n_samples={x_train.shape[0]}, n_features={x_train.shape[1]}, "
        f"n_Cs={len(c_grid)}, n_l1_ratios={len(l1_grid)}, "
        f"cv_tasks={_n_cv_tasks}, batches={_n_batches} (cores={max_cores}), "
        f"~{_est_min:.0f} min est.",
    )
    estimator = Pipeline(
        steps=[
            ("prep", preprocessor),
            (
                "model",
                LogisticRegressionCV(
                    Cs=c_grid,
                    l1_ratios=l1_grid,
                    solver="saga",
                    max_iter=20000,
                    class_weight="balanced",
                    random_state=RNG_SEED,
                    cv=inner_cv,
                    n_jobs=max_cores,
                    use_legacy_attributes=False,
                    # Keep solver output quiet; debug progress is reported at a higher level.
                    verbose=0,
                ),
            ),
        ]
    )
    # Non-convergence at extreme C values is expected during the warm-start path
    # sweep and does not affect the selected optimum — suppress the noise.
    from sklearn.exceptions import ConvergenceWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        estimator.fit(x_train, y_train)
    model = estimator.named_steps["model"]

    # Best C and l1_ratio: scalar floats with use_legacy_attributes=False.
    best_c = float(np.asarray(getattr(model, "C_", np.nan)).ravel()[0])
    l1_ratio = float(np.asarray(getattr(model, "l1_ratio_", np.nan)).ravel()[0])
    n_cs = int(np.asarray(getattr(model, "Cs_", np.array([]))).size)

    # Best CV score: scores_ is dict {class: array(n_folds, n_Cs)} (legacy) or
    # array(n_folds, n_Cs) (new API). Handle both defensively.
    best_score = float("nan")
    raw_scores = getattr(model, "scores_", None)
    if raw_scores is not None:
        try:
            if isinstance(raw_scores, dict):
                class_scores = next(iter(raw_scores.values()))
            else:
                class_scores = raw_scores
            mean_per_c = np.mean(np.asarray(class_scores, dtype=float), axis=0)
            best_score = float(np.max(mean_per_c))
        except Exception as exc:
            vlog(
                verbose,
                f"elasticnet best_score extraction failed: {exc}",
                level="debug",
            )

    vlog(
        verbose,
        f"Search complete model=elasticnet, task={task}, strategy=logisticelasticnetcv, "
        f"candidates={n_cs if n_cs > 0 else len(_C_GRID)}, best_score={best_score:.6f}",
    )

    # Replace LogisticRegressionCV with a plain LogisticRegression using the
    # selected best_c and l1_ratio. LogisticRegressionCV.fit() always re-runs
    # the full CV path sweep, so leaving it in place causes the evaluator's
    # refit call to duplicate ~19 min of work. LogisticRegression.fit() is a
    # single solver pass (~seconds).
    # Use cuML LogisticRegression for the final fit when GPU is available — a
    # single fixed-param fit avoids the precision/convergence issues that made
    # the full CV sweep prefer CPU.
    from .cuml_utils import cuml_gpu_available

    if cuml_use_gpu is not False and cuml_gpu_available():
        from cuml.linear_model import LogisticRegression as _LogisticRegression

        fixed_model = _LogisticRegression(
            C=best_c,
            l1_ratio=l1_ratio,
            penalty="elasticnet",
            class_weight="balanced",
            max_iter=20000,
            verbose=False,
        )
    else:
        from sklearn.linear_model import LogisticRegression as _LogisticRegression

        fixed_model = _LogisticRegression(
            C=best_c,
            l1_ratio=l1_ratio,
            penalty="elasticnet",
            solver="saga",
            max_iter=20000,
            class_weight="balanced",
            random_state=RNG_SEED,
        )
    estimator = Pipeline(
        steps=[("prep", estimator.named_steps["prep"]), ("model", fixed_model)]
    )

    return estimator, {
        "best_params": sanitize_best_params(
            {
                "model": fixed_model,
                "model__C": best_c,
                "model__l1_ratio": l1_ratio,
                "model__cv_folds": int(len(inner_cv)),
                "model__Cs_count": len(c_grid),
            }
        ),
        "best_score": best_score,
    }


def _fit_cuml_logistic_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    c_grid: list[float],
    l1_grid: list[float],
    verbose: bool,
    debug_grid_progress: bool = False,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit cuML LogisticRegression via GridSearchCV over a coarser C/l1_ratio grid.

    cuML LogisticRegression supports penalty="elasticnet", class_weight="balanced",
    C, and l1_ratio with the same semantics as sklearn.  No CV-path variant exists,
    so GridSearchCV drives the search (n_jobs=1 for GPU).

    DEPRECATED: not called by _fit_logistic_elasticnet_cv.  sklearn
    LogisticRegressionCV is faster and more accurate (warm-start path vs cold
    cuML float32 fits; ~30 min / 0.650 GPU vs ~20 min / 0.680 CPU at 32 cores).
    Retained for reference until a regularization-path cuML variant exists.
    """
    from .cuml_utils import get_cuml_logistic_elasticnet

    pipe = Pipeline(
        steps=[
            ("prep", preprocessor),
            ("model", get_cuml_logistic_elasticnet()),
        ]
    )
    candidates = [
        {"model__C": [float(c)], "model__l1_ratio": [float(r)]}
        for c in c_grid
        for r in l1_grid
    ]
    search = _ProgressGridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=inner_cv,
        n_jobs=1,
        refit=True,
        error_score="raise",
        progress_enabled=debug_grid_progress,
        progress_desc=f"elasticnet [{task}]",
    )
    try:
        search.fit(x_train, y_train)
    except Exception as exc:
        if _looks_like_gpu_failure(exc):
            vlog(
                verbose,
                f"cuML LogisticRegression GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            from sklearn.linear_model import (
                LogisticRegressionCV as _LogisticRegressionCV,
            )

            cpu_est = Pipeline(
                steps=[
                    ("prep", preprocessor),
                    (
                        "model",
                        _LogisticRegressionCV(
                            Cs=c_grid,
                            l1_ratios=l1_grid,
                            solver="saga",
                            max_iter=20000,
                            class_weight="balanced",
                            random_state=RNG_SEED,
                            cv=inner_cv,
                            n_jobs=1,
                            use_legacy_attributes=False,
                        ),
                    ),
                ]
            )
            cpu_est.fit(x_train, y_train)
            cpu_model = cpu_est.named_steps["model"]
            best_c = float(np.asarray(getattr(cpu_model, "C_", np.nan)).ravel()[0])
            l1_ratio = float(
                np.asarray(getattr(cpu_model, "l1_ratio_", np.nan)).ravel()[0]
            )
            best_score = float("nan")
            raw_scores = getattr(cpu_model, "scores_", None)
            if raw_scores is not None:
                try:
                    scores_arr = (
                        next(iter(raw_scores.values()))
                        if isinstance(raw_scores, dict)
                        else raw_scores
                    )
                    best_score = float(
                        np.max(np.mean(np.asarray(scores_arr, dtype=float), axis=0))
                    )
                except Exception as exc:
                    vlog(
                        verbose,
                        f"elasticnet cpu_fallback best_score extraction failed: {exc}",
                        level="debug",
                    )
            vlog(
                verbose,
                f"Search complete model=elasticnet, task={task}, strategy=logisticelasticnetcv_cpu_fallback, "
                f"candidates={len(c_grid) * len(l1_grid)}, best_score={best_score:.6f}",
            )
            return cpu_est, {
                "best_params": sanitize_best_params(
                    {
                        "model": cpu_model,
                        "model__C": best_c,
                        "model__l1_ratio": l1_ratio,
                        "model__cv_folds": int(len(inner_cv)),
                        "model__Cs_count": len(c_grid),
                    }
                ),
                "best_score": best_score,
            }
        raise

    best_params = search.best_params_
    total_candidates = len(c_grid) * len(l1_grid)
    vlog(
        verbose,
        f"Search complete model=elasticnet, task={task}, strategy=cuml_logistic_elasticnet_grid, "
        f"candidates={total_candidates}, best_score={float(search.best_score_):.6f}",
    )
    return search.best_estimator_, {
        "best_params": sanitize_best_params(best_params),
        "best_score": float(search.best_score_),
    }


def _fit_beta(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit BetaRegressor — no hyperparameters, so CV is used only for scoring.

    Uses ``sklearn.model_selection.cross_validate`` (the standard "score then
    refit" pattern) rather than a manual fold loop.

    BetaRegressor expects targets on the PSI (0, 1) scale. Logit-transformed
    targets are inverted before fitting.
    """
    y_train_beta = np.asarray(y_train, dtype=float)
    if bool(np.any(y_train_beta < 0.0) or np.any(y_train_beta > 1.0)):
        y_train_beta = _inverse_logit(y_train_beta)

    pipe = Pipeline(steps=[("prep", preprocessor), ("model", BetaRegressor())])

    # cross_validate scores the pipeline on each inner fold without leakage.
    best_score = float("nan")
    try:
        with warnings.catch_warnings():
            # Beta can fail on some folds; keep only finite fold scores.
            warnings.simplefilter("ignore", FitFailedWarning)
            cv_results = cross_validate(
                clone(pipe),
                x_train,
                y_train_beta,
                cv=inner_cv,
                scoring=scorer_name(task),
                error_score=np.nan,
                n_jobs=max_cores,
            )

        fold_scores = np.asarray(cv_results["test_score"], dtype=float)
        finite_scores = fold_scores[np.isfinite(fold_scores)]
        if finite_scores.size:
            best_score = float(np.mean(finite_scores))
        else:
            vlog(
                verbose,
                "beta inner CV produced no finite scores; continuing with full-fit beta",
            )
    except Exception as exc:
        # Inner CV is for scoring only (beta has no tuned hyperparameters).
        compact = (
            str(exc).strip().splitlines()[0]
            if str(exc).strip()
            else exc.__class__.__name__
        )
        vlog(
            verbose,
            f"beta inner CV failed ({compact}); continuing with full-fit beta",
            level="info",
        )

    # Refit on the full training fold with PSI-scale targets.
    try:
        pipe.fit(x_train, y_train_beta)
    except Exception as exc:
        compact = (
            str(exc).strip().splitlines()[0]
            if str(exc).strip()
            else exc.__class__.__name__
        )
        raise RuntimeError(f"beta full-fit failed: {compact}") from exc

    beta_model = pipe.named_steps.get("model")
    beta_risk = getattr(beta_model, "convergence_risk_", None)
    if isinstance(beta_risk, dict):
        vlog(
            verbose,
            "beta convergence precheck: "
            f"level={beta_risk.get('level')} "
            f"n={beta_risk.get('n_obs')} p={beta_risk.get('n_feat')} "
            f"p_over_n={beta_risk.get('p_over_n', float('nan')):.3f} "
            f"y_std={beta_risk.get('y_std', float('nan')):.4f} "
            f"boundary_frac={beta_risk.get('boundary_frac', float('nan')):.3f}",
            level="debug",
        )

    vlog(
        verbose,
        f"Search complete model=beta, task={task}, strategy=fixed, "
        f"candidates=1, best_score={best_score:.6f}",
    )
    return pipe, {
        "best_params": sanitize_best_params({"model": pipe.named_steps["model"]}),
        "best_score": best_score,
        "beta_convergence_risk": safe_json(beta_risk),
    }



def _slurm_available_cpu_mb() -> float | None:
    """Return available memory within the current cgroup/SLURM allocation, in MB.

    SLURM enforces limits via cgroups; psutil.virtual_memory() reads
    /proc/meminfo and therefore sees the full node RAM, not the job cap.
    We read the process's own cgroup memory controller directly.

    Returns None when not running under a memory-capped cgroup (e.g. local dev).
    """
    from pathlib import Path

    # Locate the process's own cgroup from /proc/self/cgroup.
    # cgroup v2: single line "0::/<path>"
    # cgroup v1: multiple lines "<id>:<subsystems>:<path>"; find the memory one.
    try:
        cgroup_lines = Path("/proc/self/cgroup").read_text().splitlines()
    except OSError:
        cgroup_lines = []

    # --- cgroup v2 ---
    for line in cgroup_lines:
        parts = line.split(":", 2)
        if len(parts) == 3 and parts[0] == "0":
            cg_path = parts[2].lstrip("/")
            base = Path("/sys/fs/cgroup") / cg_path
            # Walk up towards root until we find a memory.max with a real limit.
            candidate = base
            while True:
                limit_file = candidate / "memory.max"
                usage_file = candidate / "memory.current"
                try:
                    limit_txt = limit_file.read_text().strip()
                    if limit_txt not in ("max", "") and int(limit_txt) < 2**62:
                        limit = int(limit_txt)
                        usage = int(usage_file.read_text().strip())
                        return max(0.0, (limit - usage) / (1024 * 1024))
                except (OSError, ValueError):
                    pass
                parent = candidate.parent
                if parent == candidate:
                    break
                candidate = parent
            break  # v2 has only one line; stop after it

    # --- cgroup v1 ---
    for line in cgroup_lines:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        subsystems = parts[1].split(",")
        if "memory" not in subsystems:
            continue
        cg_path = parts[2].lstrip("/")
        base = Path("/sys/fs/cgroup/memory") / cg_path
        candidate = base
        while True:
            limit_file = candidate / "memory.limit_in_bytes"
            usage_file = candidate / "memory.usage_in_bytes"
            try:
                limit = int(limit_file.read_text().strip())
                if limit < 2**62:
                    usage = int(usage_file.read_text().strip())
                    return max(0.0, (limit - usage) / (1024 * 1024))
            except (OSError, ValueError):
                pass
            parent = candidate.parent
            if parent == candidate:
                break
            candidate = parent
        break

    # Fallback: SLURM_MEM_PER_NODE is in MB; subtract a small buffer for OS overhead.
    mem_str = os.environ.get("SLURM_MEM_PER_NODE")
    if mem_str:
        try:
            return max(0.0, float(mem_str) - 2048)  # keep 2 GB for OS
        except ValueError:
            pass

    return None


def _patch_psutil_for_slurm() -> None:
    """Permanently patch psutil.virtual_memory to respect the cgroup/SLURM memory cap.

    TabICL's InferenceManager.get_available_cpu_memory() calls
    psutil.virtual_memory().available during predict(), which reads /proc/meminfo
    and sees the full node RAM — not the SLURM job's cgroup limit.  A context
    manager cannot fix this because predict() is called from a different call
    stack (evaluate_outer_fold) long after _fit_tabicl() returns.

    Applying the patch once here is safe: the patched function re-reads the
    cgroup on every call so it stays accurate as memory usage changes, and the
    cap is always min(real, cgroup_limit) so it can never overstate available RAM.
    No-op when not under a memory-capped cgroup.
    """
    import psutil

    avail_mb = _slurm_available_cpu_mb()
    if avail_mb is None:
        return  # not under a cgroup limit; nothing to do

    if getattr(psutil.virtual_memory, "_slurm_patched", False):
        return  # already patched; don't double-wrap

    _real_vmem = psutil.virtual_memory
    cap_bytes = int(avail_mb * 1024 * 1024)

    def _patched_vmem():
        real = _real_vmem()
        # Re-read cgroup available on every call so the cap tracks actual usage.
        current_cap = _slurm_available_cpu_mb()
        if current_cap is not None:
            cap = int(current_cap * 1024 * 1024)
        else:
            cap = cap_bytes
        return real._replace(available=min(real.available, cap))

    _patched_vmem._slurm_patched = True  # type: ignore[attr-defined]
    psutil.virtual_memory = _patched_vmem


def _fit_tabicl(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    verbose: bool,
    tabicl_n_estimators: int = 1,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit TabICLv2 directly on the outer training fold — no inner CV needed.

    TabICL is a pre-trained in-context learning model with no hyperparameters
    to tune. fit() stores the training data; predict() runs a transformer
    forward pass using that data as context. Checkpoints are auto-downloaded
    from HuggingFace on first use.
    """
    from tabicl import TabICLClassifier, TabICLRegressor

    # Patch psutil permanently so the SLURM cgroup limit is visible to TabICL
    # during both fit() and predict() (predict() is called later from evaluate_outer_fold).
    _patch_psutil_for_slurm()

    common_kwargs: dict[str, Any] = dict(
        n_estimators=tabicl_n_estimators,
        verbose=verbose,
        random_state=RNG_SEED,
        batch_size=1,
        offload_mode="auto",  # psutil is patched above; auto decision is now cgroup-aware
        disk_offload_dir=os.path.join(
            os.environ.get("TMPDIR", "/localscratch"), "tabicl_offload"
        ),  # absolute path — relative paths silently fail under SLURM, disabling disk fallback
    )

    if task == "regression":
        model = TabICLRegressor(
            checkpoint_version="tabicl-regressor-v2-20260212.ckpt",
            **common_kwargs,
        )
    else:
        model = TabICLClassifier(
            checkpoint_version="tabicl-classifier-v2-20260212.ckpt",
            **common_kwargs,
        )

    pipe = Pipeline(steps=[("prep", preprocessor), ("model", model)])

    try:
        pipe.fit(x_train, y_train)
    except Exception as exc:
        compact = (
            str(exc).strip().splitlines()[0]
            if str(exc).strip()
            else exc.__class__.__name__
        )
        raise RuntimeError(f"tabicl fit failed: {compact}") from exc

    vlog(
        verbose,
        f"Search complete model=tabicl, task={task}, strategy=fixed, candidates=1",
    )
    return pipe, {
        "best_params": sanitize_best_params({"model": pipe.named_steps["model"]}),
        "best_score": float("nan"),
    }


def _looks_like_gpu_failure(exc: Exception) -> bool:
    """Return True when an exception message suggests a CUDA/GPU problem."""
    msg = str(exc).lower()
    return any(
        token in msg
        for token in ("cuda", "gpu", "out of memory", "oom", "device", "cublas", "nccl")
    )


def _compute_es_splits(
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Compute per-fold early-stopping (ES) train/val index pairs.

    For each inner fold i, the ES val is the smallest scoring-val from any
    *other* inner fold — it is already inside fold i's train set (by
    construction) and is disjoint from fold i's scoring val, so no leakage
    between the ES signal and the CV score.  The ES train is fold i's train
    minus the borrowed ES val, keeping training data as large as possible.

    Returns a list of (es_train_idx, es_val_idx) arrays, one per fold.
    """
    n = len(inner_cv)
    result: list[tuple[np.ndarray, np.ndarray]] = []
    for i, (train_idx, _) in enumerate(inner_cv):
        es_val_idx = min(
            (inner_cv[j][1] for j in range(n) if j != i),
            key=len,
        )
        es_val_set = set(es_val_idx.tolist())
        es_train_idx = np.array(
            [idx for idx in train_idx if idx not in es_val_set], dtype=int
        )
        result.append((es_train_idx, es_val_idx))
    return result


class _ContextAwareCV:
    """CV splitter that injects a clean ES val before each fold.

    For mlp / xgb, computes per-fold early-stopping splits via
    _compute_es_splits:
      - ES val: the smallest scoring-val from any *other* inner fold
        (already inside the current fold's train; disjoint from the
        scoring val, so no leakage between ES signal and CV score).
      - ES train: current fold's train minus the borrowed ES val.

    Stores (X_es_val_raw, y_es_val) in _es_raw_val_context so
    _ESPipeline can preprocess and pass to the model for early stopping.
    Yields (es_train_idx, scoring_val_idx) to sklearn so the scoring val
    is never touched by early stopping.
    """

    def __init__(
        self,
        inner_cv: list[tuple[np.ndarray, np.ndarray]],
        model_name: str,
    ) -> None:
        """Initialize a _ContextAwareCV instance."""
        self._inner_cv = inner_cv
        self._model_name = model_name
        self._es_splits: list[tuple[np.ndarray, np.ndarray]] | None = None
        if model_name in {"mlp", "xgb", "lgbm"} and len(inner_cv) >= 2:
            self._es_splits = _compute_es_splits(inner_cv)

    def split(self, X: Any, y: Any = None, groups: Any = None):  # type: ignore[override]
        """Yield (es_train_idx, scoring_val_idx) for each fold."""
        for i, (train_idx, val_idx) in enumerate(self._inner_cv):
            if self._model_name in {"mlp", "xgb", "lgbm"} and self._es_splits is not None:
                es_train_idx, es_val_idx = self._es_splits[i]
                _es_raw_val_context.X_val = (
                    X.iloc[es_val_idx] if hasattr(X, "iloc") else X[es_val_idx]
                )
                _es_raw_val_context.y_val = y[es_val_idx] if y is not None else None
                yield es_train_idx, val_idx
            else:
                yield train_idx, val_idx
        # Clear after all folds to avoid stale state.
        _es_raw_val_context.X_val = None
        _es_raw_val_context.y_val = None

    def get_n_splits(self, X: Any = None, y: Any = None, groups: Any = None) -> int:
        """Get n splits."""
        return len(self._inner_cv)


_TREE_ES_ROUNDS = (
    30  # patience for per-trial XGB early stopping (consecutive non-improving rounds)
)
_XGB_ES_N_MAX = 2000  # max n_estimators for XGB post-search probe
_XGB_ES_N_MIN = 50  # floor on XGB optimal n_estimators (degenerate below this)
_LGBM_ES_N_MAX = 2000  # max n_estimators for LightGBM post-search probe (mirrors XGB's; leaf-wise vs depth-wise tree counts aren't strictly comparable, may need retuning)
_LGBM_ES_N_MIN = 50  # floor on LightGBM optimal n_estimators (degenerate below this)
_MLP_ES_N_MAX = 500  # max epochs for MLP post-search probe (hard ceiling)
# MLP_ES_PATIENCE imported from grids (defined alongside _mlp_base_model so one source of truth)


def _lgbm_classification_sample_weight(y: np.ndarray) -> np.ndarray:
    """Per-sample class-imbalance weight for LightGBM classification fits.

    LightGBM's CUDA backend has a broken internal class-reweighting code path
    -- both `scale_pos_weight=` and `is_unbalance=True` corrupt gradients on
    `device_type="cuda"`, collapsing training to a root-only tree (confirmed
    2026-07-14, see CLAUDE.md). Passing the equivalent weight as per-sample
    `sample_weight=` at fit time sidesteps it (computed in Python, outside
    CUDA's broken reweighting path).
    """
    n_pos = int(np.sum(y == 1))
    n_neg = int(len(y) - n_pos)
    ratio = float(n_neg) / max(1, n_pos)
    return np.where(y == 1, ratio, 1.0).astype(np.float64)


class _ESPipeline(Pipeline):
    """Pipeline subclass that enables per-fold early stopping for XGB/LightGBM/MLP.

    When _es_raw_val_context holds raw val data (set by _ContextAwareCV),
    _ESPipeline preprocesses it with the just-fitted preprocessor and:
    - For XGB: passes eval_set + early_stopping_rounds to model.fit()
    - For LightGBM: passes eval_set + a lightgbm.early_stopping callback to
      model.fit() (LightGBM 4.x removed the early_stopping_rounds/verbose fit
      kwargs in favor of callbacks=[...])
    - For MLP: stores preprocessed val in _es_pp_val_context for the model to read

    Falls back to standard Pipeline.fit() when no val data is available
    (e.g. during OptunaSearchCV's final refit on all x_train).
    """

    def __init__(self, steps, *, model_name: str = "", **kwargs):
        """Initialize a _ESPipeline instance."""
        super().__init__(steps, **kwargs)
        self.model_name = model_name
        # Expose estimator type from final estimator
        final_est = self._final_estimator
        if hasattr(final_est, "_estimator_type"):
            self._estimator_type = final_est._estimator_type

    def fit(self, X, y=None, **params):
        # Fit all intermediate steps (preprocessor) via sklearn's internal _fit().
        # _fit() fits each step up to (not including) the last and returns Xt.
        """Fit the model using the provided training data."""
        routed_params = self._check_method_params(method="fit", props=params)
        Xt = self._fit(X, y, routed_params, raw_params=params)

        # Preprocess val data using the now-fitted preprocessor(s).
        X_val_raw = getattr(_es_raw_val_context, "X_val", None)
        y_val = getattr(_es_raw_val_context, "y_val", None)
        X_val_pp = None
        if X_val_raw is not None:
            X_val_pp = X_val_raw
            for _, _, transformer in self._iter(
                with_final=False, filter_passthrough=False
            ):
                X_val_pp = transformer.transform(X_val_pp)

        # Fit the final estimator with early stopping.
        model = self._final_estimator
        last_name = self.steps[-1][0]
        last_params = self._get_metadata_for_step(
            step_idx=len(self) - 1,
            step_params=routed_params[last_name],
            all_params=params,
        )
        fit_params = last_params.get("fit", {})

        if X_val_pp is not None and y_val is not None and self.model_name == "xgb":
            try:
                # early_stopping_rounds is a constructor param in XGBoost >= 2.x.
                model.set_params(early_stopping_rounds=_TREE_ES_ROUNDS)
                model.fit(
                    Xt,
                    y,
                    eval_set=[(X_val_pp, y_val)],
                    verbose=False,
                    **fit_params,
                )
            except Exception as exc:
                vlog(
                    True,
                    f"XGB early stopping failed, falling back to standard fit: {exc}",
                    level="debug",
                )
                model.fit(Xt, y, **fit_params)
        elif X_val_pp is not None and y_val is not None and self.model_name == "lgbm":
            lgbm_fit_params = dict(fit_params)
            if model.__class__.__name__ == "LGBMClassifier":
                lgbm_fit_params["sample_weight"] = _lgbm_classification_sample_weight(y)
            try:
                import lightgbm

                model.fit(
                    Xt,
                    y,
                    eval_set=[(X_val_pp, y_val)],
                    callbacks=[
                        lightgbm.early_stopping(
                            stopping_rounds=_TREE_ES_ROUNDS, verbose=False
                        )
                    ],
                    **lgbm_fit_params,
                )
            except Exception as exc:
                vlog(
                    True,
                    f"LightGBM early stopping failed, falling back to standard fit: {exc}",
                    level="debug",
                )
                model.fit(Xt, y, **lgbm_fit_params)
        elif X_val_pp is not None and y_val is not None and self.model_name == "mlp":
            _es_pp_val_context.X_val_pp = X_val_pp
            _es_pp_val_context.y_val = y_val
            try:
                model.fit(Xt, y, **fit_params)
            finally:
                _es_pp_val_context.X_val_pp = None
                _es_pp_val_context.y_val = None
        elif self.model_name == "lgbm" and model.__class__.__name__ == "LGBMClassifier":
            model.fit(
                Xt, y,
                sample_weight=_lgbm_classification_sample_weight(y),
                **fit_params,
            )
        else:
            model.fit(Xt, y, **fit_params)

        return self


def _build_grid_search(
    pipe: Pipeline,
    model_name: str,
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    search_n_jobs: int,
    param_grid_size: int,
    search_strategy: str,
    lhs_scale_mode: str,
    xgb_use_gpu: bool | None,
    cuml_use_gpu: bool | None,
    verbose: bool,
    model_n_jobs: int = 1,
    debug_grid_progress: bool = False,
    lgbm_use_gpu: bool | None = None,
    random_seed: int = RNG_SEED,
) -> GridSearchCV:
    """Construct a GridSearchCV object with the appropriate candidate list."""
    # NOTE: for lgbm, scale_pos_weight is computed here for backward-compatible
    # threading through build_param_candidates, but _build_lgbm_base_model no
    # longer forwards it to the constructor (unsafe on device_type="cuda" --
    # see _lgbm_classification_sample_weight). Imbalance is applied via
    # sample_weight at fit time instead. Still used normally for "rf" (XGBRF,
    # unaffected).
    scale_pos_weight = None
    if model_name in {"rf", "lgbm"} and task == "classification":
        n_pos = int(np.sum(y_train == 1))
        n_neg = int(np.sum(y_train == 0))
        scale_pos_weight = float(n_neg) / max(1, n_pos)

    candidates = build_param_candidates(
        model_name=model_name,
        task=task,
        n_samples=x_train.shape[0],
        n_features=x_train.shape[1],
        budget=param_grid_size,
        strategy=search_strategy,
        lhs_scale_mode=lhs_scale_mode,
        xgb_use_gpu=xgb_use_gpu,
        cuml_use_gpu=cuml_use_gpu,
        lgbm_use_gpu=lgbm_use_gpu,
        verbose=verbose,
        model_n_jobs=model_n_jobs,
        scale_pos_weight=scale_pos_weight,
        seed=random_seed,
    )

    # For XGB with early stopping, n_estimators is determined by ES — not by the
    # search.  Fix it to _XGB_ES_N_MAX on every candidate's base model and remove
    # it from the grid so the search only tunes the remaining hyperparameters.
    if model_name == "xgb":
        for c in candidates:
            c.pop("model__n_estimators", None)
            for m in c.get("model", []):
                m.set_params(n_estimators=_XGB_ES_N_MAX)
        vlog(
            verbose,
            f"XGB grid search: fixed n_estimators={_XGB_ES_N_MAX} (determined by ES); "
            f"removed from search grid",
            level="debug",
        )

    # Same treatment for LightGBM: n_estimators comes from post-search ES.
    if model_name == "lgbm":
        for c in candidates:
            c.pop("model__n_estimators", None)
            for m in c.get("model", []):
                m.set_params(n_estimators=_LGBM_ES_N_MAX)
        vlog(
            verbose,
            f"LGBM grid search: fixed n_estimators={_LGBM_ES_N_MAX} (determined by ES); "
            f"removed from search grid",
            level="debug",
        )

    # Wrap the pipeline with _ESPipeline for models that need early stopping.
    if model_name in {"xgb", "mlp", "lgbm"}:
        pipe = _ESPipeline(pipe.steps, model_name=model_name)

    # Activate _ContextAwareCV for models that need inner-fold validation data.
    cv_splitter = inner_cv
    if model_name in {"xgb", "mlp", "lgbm"}:
        cv_splitter = _ContextAwareCV(inner_cv, model_name)

    return _ProgressGridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=cv_splitter,
        n_jobs=search_n_jobs,
        pre_dispatch=search_n_jobs,
        refit=True,
        error_score="raise",
        progress_enabled=debug_grid_progress,
        progress_desc=f"grid {model_name} [{task}]",
    )


def _fit_tree_early_stopping(
    model_name: str,
    task: str,
    fitted_pipe: Pipeline,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    verbose: bool,
    reserve_calibration_holdout: bool = False,
) -> tuple[Pipeline, bool, int | None, tuple[pd.DataFrame, np.ndarray] | None]:
    """Post-search refit with early stopping on final model (XGB, LightGBM, MLP only).

    Uses the globally smallest inner-fold val as the ES hold-out, then trains
    the probe model on ALL outer training data except that hold-out:
      - ES val: the inner-fold val with fewest samples (smallest group).
      - Probe train: all outer train indices minus ES val.
    This maximises the probe training set so it best approximates the final fit,
    while keeping the ES signal on a partition that was never the scoring val
    for the same fold it was trained on.

    Only applies early stopping to XGBoost, LightGBM, and MLP (boosted/iterative
    models). RandomForest (non-boosted) doesn't support early stopping with
    eval_set, so it just uses standard fit with the n_estimators from the search.

    For XGB: Uses eval_set with early_stopping_rounds.
    For LightGBM: Uses eval_set with a lightgbm.early_stopping callback (4.x
    removed the early_stopping_rounds/verbose fit kwargs).
    For MLP: Stores validation data in _es_pp_val_context for the model to read.
    For RF: Standard fit (no early stopping available for random forests).

    reserve_calibration_holdout (XGB and LightGBM only): instead of refitting
    the final model on ALL outer-train data, keep the ES-val partition held
    out and return it (raw, pre-preprocessing) as the 4th tuple element so the
    caller can fit ``CalibratedClassifierCV(..., estimator=FrozenEstimator(...))``
    on a slice that was never seen by the fitted model — unlike calibrating via
    ``cv=inner_splits``, which reuses the exact folds that already selected the
    hyperparameters. Costs the final model that one held-out slice of training
    data; ignored for model_name not in {"xgb", "lgbm"} (no other model routes
    through this reservation path today).

    Returns (pipe, applied, optimal_n, calibration_holdout).
    """
    try:
        prep_fitted = fitted_pipe.named_steps["prep"]
        best_model = fitted_pipe.named_steps["model"]

        # Use the globally smallest inner-fold val as ES hold-out, then train the
        # probe on ALL outer training data except that hold-out.  This maximises
        # the probe training set (best approximation of the final fit) while
        # keeping the ES signal on a clean, separate partition.
        es_val_idx = min((val_idx for _, val_idx in inner_cv), key=len)
        es_val_set = set(es_val_idx.tolist())
        all_train_idx = np.arange(len(x_train))
        probe_tr_idx = all_train_idx[
            np.array([i not in es_val_set for i in all_train_idx])
        ]
        X_final_train_t = prep_fitted.transform(x_train.iloc[probe_tr_idx])
        X_final_val_t = prep_fitted.transform(x_train.iloc[es_val_idx])
        y_final_train = y_train[probe_tr_idx]
        y_final_val = y_train[es_val_idx]

        model_final = clone(best_model)
        optimal_n = None
        # Populated in the lgbm branch below (from y_final_train); reused by the
        # calibration-holdout refit further down, which fits on the same
        # y_final_train. The normal (non-calibration) final refit uses the full
        # y_train instead and gets its own freshly-computed weight array there.
        lgbm_sw: np.ndarray | None = None

        if model_name == "xgb":
            # XGBoost (boosted): try early stopping, fall back to standard fit if unsupported
            model_final.set_params(n_estimators=_XGB_ES_N_MAX)
            if y_final_val is not None:
                vlog(
                    verbose,
                    f"XGB post-search ES: probe_train_n={len(X_final_train_t)}, es_val_n={len(X_final_val_t)}",
                    level="debug",
                )
                try:
                    # early_stopping_rounds is a constructor param in XGBoost >= 2.x.
                    model_final.set_params(early_stopping_rounds=_TREE_ES_ROUNDS)
                    model_final.fit(
                        X_final_train_t,
                        y_final_train,
                        eval_set=[(X_final_val_t, y_final_val)],
                        verbose=False,
                    )
                    optimal_n = max(_XGB_ES_N_MIN, int(model_final.best_iteration) + 1)
                except Exception as exc:
                    vlog(
                        verbose,
                        f"XGB post-search early stopping failed, falling back: {exc}",
                        level="debug",
                    )
                    model_final.fit(X_final_train_t, y_final_train)
                    optimal_n = _XGB_ES_N_MAX
            else:
                model_final.fit(X_final_train_t, y_final_train)
                optimal_n = _XGB_ES_N_MAX
        elif model_name == "lgbm":
            # LightGBM (boosted): try early stopping, fall back to standard fit if unsupported
            model_final.set_params(n_estimators=_LGBM_ES_N_MAX)
            # scale_pos_weight/is_unbalance are unsafe on device_type="cuda" (see
            # _lgbm_classification_sample_weight's docstring) -- classification
            # imbalance is applied via sample_weight at every fit call below instead.
            lgbm_sw = (
                _lgbm_classification_sample_weight(y_final_train)
                if task == "classification"
                else None
            )
            if y_final_val is not None:
                vlog(
                    verbose,
                    f"LGBM post-search ES: probe_train_n={len(X_final_train_t)}, es_val_n={len(X_final_val_t)}",
                    level="debug",
                )
                try:
                    import lightgbm

                    model_final.fit(
                        X_final_train_t,
                        y_final_train,
                        eval_set=[(X_final_val_t, y_final_val)],
                        callbacks=[
                            lightgbm.early_stopping(
                                stopping_rounds=_TREE_ES_ROUNDS, verbose=False
                            )
                        ],
                        sample_weight=lgbm_sw,
                    )
                    # LightGBM's sklearn wrapper exposes best_iteration_ (trailing
                    # underscore) -- unlike XGBoost's best_iteration (no underscore).
                    optimal_n = max(
                        _LGBM_ES_N_MIN, int(model_final.best_iteration_) + 1
                    )
                except Exception as exc:
                    vlog(
                        verbose,
                        f"LGBM post-search early stopping failed, falling back: {exc}",
                        level="debug",
                    )
                    model_final.fit(X_final_train_t, y_final_train, sample_weight=lgbm_sw)
                    optimal_n = _LGBM_ES_N_MAX
            else:
                model_final.fit(X_final_train_t, y_final_train, sample_weight=lgbm_sw)
                optimal_n = _LGBM_ES_N_MAX
        elif model_name == "mlp":
            # MLP: inject val context so early stopping uses a proper group partition,
            # not the sequential val_fraction fallback (which is typically OOD).
            if y_final_val is not None:
                _es_pp_val_context.X_val_pp = X_final_val_t
                _es_pp_val_context.y_val = y_final_val
                vlog(
                    verbose,
                    f"MLP post-search ES: probe_train_n={len(X_final_train_t)}, es_val_n={len(X_final_val_t)}",
                    level="debug",
                )
                try:
                    model_final.fit(X_final_train_t, y_final_train)
                finally:
                    _es_pp_val_context.X_val_pp = None
                    _es_pp_val_context.y_val = None
            else:
                model_final.fit(X_final_train_t, y_final_train)
            # Use discovered epoch count for final refit; cap at _MLP_ES_N_MAX.
            # No floor: if ES genuinely converged early, honour that.
            optimal_n = min(
                _MLP_ES_N_MAX,
                int(getattr(model_final, "n_epochs_trained_", _MLP_ES_N_MAX)),
            )
        else:
            # Fallback for unknown model types
            model_final.fit(X_final_train_t, y_final_train)
            optimal_n = _MLP_ES_N_MAX

        # Refit final model on all training data using the discovered optimal depth.
        # For XGB: set n_estimators=optimal_n and clear early_stopping_rounds so the
        #   fit runs for exactly optimal_n trees without needing an eval_set.
        # For MLP: set max_epochs=optimal_n and disable weight restoration — the probe
        #   already identified the right epoch count, and there is no clean val set for
        #   all-data training (es_val is a subset of x_train, so it would be in-sample).
        if optimal_n is not None and model_name == "xgb":
            model_final.set_params(n_estimators=optimal_n, early_stopping_rounds=None)
        elif optimal_n is not None and model_name == "lgbm":
            # No early_stopping_rounds constructor param to clear (unlike XGB) --
            # LightGBM's ES was applied purely via a fit-time callback.
            model_final.set_params(n_estimators=optimal_n)
        elif optimal_n is not None and model_name == "mlp":
            model_final.set_params(max_epochs=optimal_n, restore_best_weights=False)

        calibration_holdout: tuple[pd.DataFrame, np.ndarray] | None = None
        if reserve_calibration_holdout and model_name in {"xgb", "lgbm"}:
            # Keep the ES-val slice held out instead of folding it back in, so
            # it can serve as a leakage-free cv="prefit" calibration set below.
            if model_name == "lgbm" and task == "classification":
                model_final.fit(X_final_train_t, y_final_train, sample_weight=lgbm_sw)
            else:
                model_final.fit(X_final_train_t, y_final_train)
            calibration_holdout = (x_train.iloc[es_val_idx], y_final_val)
        else:
            X_all_train_t = prep_fitted.transform(x_train)
            if model_name == "lgbm" and task == "classification":
                model_final.fit(
                    X_all_train_t, y_train,
                    sample_weight=_lgbm_classification_sample_weight(y_train),
                )
            else:
                model_final.fit(X_all_train_t, y_train)

        pipe_final = Pipeline([("prep", prep_fitted), ("model", model_final)])

        vlog(
            verbose,
            f"{model_name.upper()} post-search ES: optimal_n={optimal_n} "
            f"(final-fold validation)"
            + (
                ", calibration_holdout reserved"
                if calibration_holdout is not None
                else ""
            ),
        )
        return pipe_final, True, optimal_n, calibration_holdout

    except Exception as exc:
        vlog(
            verbose, f"{model_name.upper()} post-search ES failed ({exc})", level="info"
        )
        return fitted_pipe, False, None, None


def _log_early_stopping_fold_sizes(
    model_name: str,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    *,
    verbose: bool,
) -> None:
    """Log a compact summary of early-stopping fold sizes once per search.

    Shows three sizes per fold:
      es_train_n   – samples the model trains on (fold train minus ES val)
      scoring_val_n – held-out samples used for CV scoring (untouched by ES)
      es_val_n     – borrowed hold-out used only for early stopping
    """
    if model_name not in {"xgb", "mlp", "lgbm"} or not inner_cv:
        return

    es_splits = _compute_es_splits(inner_cv)
    fold_summaries = "; ".join(
        f"split{i} es_train_n={len(es_tr)}, scoring_val_n={len(inner_cv[i - 1][1])}, es_val_n={len(es_va)}"
        for i, (es_tr, es_va) in enumerate(es_splits, start=1)
    )
    vlog(
        verbose,
        f"{model_name.upper()} early stopping fold sizes: {fold_summaries}",
        level="debug",
    )


def _fit_grid_search(
    task: str,
    model_name: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    pipe: Pipeline,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    param_grid_size: int,
    search_strategy: str,
    lhs_scale_mode: str,
    xgb_use_gpu: bool | None,
    cuml_use_gpu: bool | None,
    verbose: bool,
    debug_grid_progress: bool,
    calibrate: bool = False,
    lgbm_use_gpu: bool | None = None,
    random_seed: int = RNG_SEED,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit all non-ElasticNet, non-Beta models via GridSearchCV.

    Handles:
    - Parallelism: limits XGBoost/LightGBM/cuML GPU jobs to 1 to avoid VRAM contention.
    - GPU retry: on CUDA failure, automatically falls back to CPU.
    - XGBoost/LightGBM post-search refit: after tuning with a reduced
      n_estimators budget, the final estimator is re-trained with a larger
      tree budget.
    - calibrate=True (XGB and LightGBM only): reserves the ES-val slice as a
      dedicated calibration holdout instead of folding it into the final refit.
    """
    # Limit concurrent GPU jobs to avoid VRAM contention.
    # For RF, split cores between search-level and tree-level parallelism:
    # model_n_jobs gets sqrt(max_cores) cores per forest; search_n_jobs gets the rest.
    if model_name == "mlp":
        # Detect any available accelerator to decide parallelism budget.
        _mlp_has_gpu = False
        try:
            import torch as _torch

            _mlp_has_gpu = _torch.cuda.is_available()
            try:
                _mlp_has_gpu = _mlp_has_gpu or _torch.accelerator.is_available()
            except AttributeError:
                vlog(
                    verbose,
                    "torch.accelerator not available (older torch); using cuda/mps only",
                    level="debug",
                )
            if hasattr(_torch.backends, "mps"):
                _mlp_has_gpu = _mlp_has_gpu or _torch.backends.mps.is_available()
        except ImportError:
            vlog(verbose, "torch not installed; MLP will run on CPU", level="debug")
        # Single search job when a GPU is available to prevent VRAM contention
        # across parallel GridSearchCV workers.
        search_n_jobs = 1 if _mlp_has_gpu else max_cores
        model_n_jobs = 1
    elif (
        (model_name in {"xgb", "rf"} and bool(xgb_use_gpu))
        or (model_name == "svm" and bool(cuml_use_gpu))
        or (model_name == "lgbm" and bool(lgbm_use_gpu))
    ):
        search_n_jobs = 1
        model_n_jobs = 1
    elif model_name in {"rf", "xgb", "lgbm"} and max_cores > 1:
        # Split budget between search-level and model-level parallelism.
        # RF/LightGBM use joblib/native tree-building threads; XGB CPU uses
        # its own thread pool.
        model_n_jobs = max(1, int(max_cores**0.5))
        search_n_jobs = max(1, max_cores // model_n_jobs)
    else:
        search_n_jobs = max_cores
        model_n_jobs = 1

    search = _build_grid_search(
        pipe=pipe,
        model_name=model_name,
        task=task,
        x_train=x_train,
        y_train=y_train,
        inner_cv=inner_cv,
        search_n_jobs=search_n_jobs,
        param_grid_size=param_grid_size,
        search_strategy=search_strategy,
        lhs_scale_mode=lhs_scale_mode,
        xgb_use_gpu=xgb_use_gpu,
        cuml_use_gpu=cuml_use_gpu,
        lgbm_use_gpu=lgbm_use_gpu,
        verbose=verbose,
        model_n_jobs=model_n_jobs,
        debug_grid_progress=debug_grid_progress,
        random_seed=random_seed,
    )

    def _fit_search(search_obj: GridSearchCV) -> None:
        """Internal helper for fit search."""
        if search_n_jobs > 1:
            # Thread backend avoids occasional loky worker-stop warnings.
            with parallel_backend("threading", n_jobs=search_n_jobs):
                search_obj.fit(x_train, y_train)
            return
        search_obj.fit(x_train, y_train)

    _log_early_stopping_fold_sizes(model_name, inner_cv, verbose=verbose)
    try:
        _fit_search(search)
    except Exception as exc:
        # On GPU failure, rebuild the search with CPU and retry once.
        if (
            model_name in {"xgb", "rf"}
            and bool(xgb_use_gpu)
            and _looks_like_gpu_failure(exc)
        ):
            vlog(
                verbose,
                f"{model_name.upper()} GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            search = _build_grid_search(
                pipe=pipe,
                model_name=model_name,
                task=task,
                x_train=x_train,
                y_train=y_train,
                inner_cv=inner_cv,
                search_n_jobs=max_cores,
                param_grid_size=param_grid_size,
                search_strategy=search_strategy,
                lhs_scale_mode=lhs_scale_mode,
                xgb_use_gpu=False,  # force CPU
                cuml_use_gpu=cuml_use_gpu,
                lgbm_use_gpu=lgbm_use_gpu,
                verbose=verbose,
                debug_grid_progress=debug_grid_progress,
                random_seed=random_seed,
            )
            _fit_search(search)
        elif (
            model_name == "svm" and bool(cuml_use_gpu) and _looks_like_gpu_failure(exc)
        ):
            vlog(
                verbose,
                f"cuML SVM GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            search = _build_grid_search(
                pipe=pipe,
                model_name=model_name,
                task=task,
                x_train=x_train,
                y_train=y_train,
                inner_cv=inner_cv,
                search_n_jobs=max_cores,
                param_grid_size=param_grid_size,
                search_strategy=search_strategy,
                lhs_scale_mode=lhs_scale_mode,
                xgb_use_gpu=xgb_use_gpu,
                cuml_use_gpu=False,  # force sklearn CPU path
                lgbm_use_gpu=lgbm_use_gpu,
                verbose=verbose,
                debug_grid_progress=debug_grid_progress,
                random_seed=random_seed,
            )
            _fit_search(search)
        elif (
            model_name == "lgbm" and bool(lgbm_use_gpu) and _looks_like_gpu_failure(exc)
        ):
            vlog(
                verbose,
                f"LGBM GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            search = _build_grid_search(
                pipe=pipe,
                model_name=model_name,
                task=task,
                x_train=x_train,
                y_train=y_train,
                inner_cv=inner_cv,
                search_n_jobs=max_cores,
                param_grid_size=param_grid_size,
                search_strategy=search_strategy,
                lhs_scale_mode=lhs_scale_mode,
                xgb_use_gpu=xgb_use_gpu,
                cuml_use_gpu=cuml_use_gpu,
                lgbm_use_gpu=False,  # force CPU
                verbose=verbose,
                debug_grid_progress=debug_grid_progress,
                random_seed=random_seed,
            )
            _fit_search(search)
        else:
            raise

    # XGB/LightGBM/MLP: post-search refit with early stopping to find optimal n_estimators/epochs.
    tree_es_applied = False
    tree_es_n_estimators: int | None = None
    calibration_holdout: tuple[pd.DataFrame, np.ndarray] | None = None
    if model_name in {"xgb", "mlp", "lgbm"}:
        pipe_es, tree_es_applied, tree_es_n_estimators, calibration_holdout = (
            _fit_tree_early_stopping(
                model_name=model_name,
                task=task,
                fitted_pipe=search.best_estimator_,
                x_train=x_train,
                y_train=y_train,
                inner_cv=inner_cv,
                verbose=verbose,
                reserve_calibration_holdout=calibrate and model_name in {"xgb", "lgbm"},
            )
        )
        search.best_estimator_ = pipe_es

    # Compute total candidate count for logging.
    use_lhs = _use_lhs_strategy(search_strategy, model_name)
    grid = search.param_grid  # list of dicts
    if use_lhs:
        total_candidates = max(1, int(param_grid_size))
    else:
        total_candidates = sum(
            math.prod(len(v) for k, v in g.items() if k != "model")
            for g in (grid if isinstance(grid, list) else [grid])
        )

    vlog(
        verbose,
        f"Search complete model={model_name}, task={task}, "
        f"strategy={'lhs' if use_lhs else 'grid'}, "
        f"candidates={total_candidates}, best_score={float(search.best_score_):.6f}",
    )

    # Switch best estimator to CPU to avoid device-mismatch warnings downstream.
    _set_xgb_cpu_predictor_for_inference(search.best_estimator_)
    _set_lgbm_cpu_predictor_for_inference(search.best_estimator_)

    tuning_info: dict[str, Any] = {
        "best_params": sanitize_best_params(search.best_params_),
        "best_score": float(search.best_score_),
    }
    if model_name in {"xgb", "mlp", "rf", "lgbm"}:
        tuning_info["tree_es_applied"] = tree_es_applied
        tuning_info["tree_es_n_estimators"] = tree_es_n_estimators
    if calibration_holdout is not None:
        tuning_info["calibration_holdout"] = calibration_holdout

    try:
        cv_res = search.cv_results_
        split_keys = sorted(
            k for k in cv_res if k.startswith("split") and k.endswith("_test_score")
        )
        n = len(cv_res["params"])
        tuning_info["cv_results"] = {
            "n_candidates": n,
            "mean_test_score": [
                float(v) if np.isfinite(float(v)) else None
                for v in cv_res["mean_test_score"]
            ],
            "std_test_score": [
                float(v) if np.isfinite(float(v)) else None
                for v in cv_res["std_test_score"]
            ],
            "rank_test_score": [int(v) for v in cv_res["rank_test_score"]],
            "params_json": [
                json.dumps(sanitize_best_params(p), sort_keys=True, default=str)
                for p in cv_res["params"]
            ],
            "split_keys": split_keys,
            "split_scores": {
                k: [float(v) if np.isfinite(float(v)) else None for v in cv_res[k]]
                for k in split_keys
            },
        }
    except Exception as exc:
        vlog(verbose, f"grid search cv_results extraction failed: {exc}", level="debug")

    return search.best_estimator_, tuning_info


def _fit_optuna_search(
    task: str,
    model_name: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    pipe: Pipeline,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    param_grid_size: int,
    optuna_sampler: str,
    optuna_n_startup_trials: int,
    optuna_multivariate: bool,
    xgb_use_gpu: bool | None,
    cuml_use_gpu: bool | None,
    verbose: bool,
    optuna_wandb_callback: bool = False,
    calibrate: bool = False,
    lgbm_use_gpu: bool | None = None,
    random_seed: int = RNG_SEED,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit a model via Optuna Bayesian optimisation (TPE or CMA-ES).

    Uses ``OptunaSearchCV`` as a drop-in for GridSearchCV.  GPU retry logic
    and XGBoost/LightGBM post-search refit are preserved from ``_fit_grid_search``.
    calibrate=True (XGB and LightGBM only): reserves the ES-val slice as a
    dedicated calibration holdout instead of folding it into the final refit.
    """
    try:
        import optuna
        from optuna_integration.sklearn import OptunaSearchCV
    except ImportError as exc:
        raise ImportError(
            "optuna_backend=True requires optuna and optuna-integration[sklearn]. "
            "Install with: pip install 'optuna-integration[sklearn]>=3.4'"
        ) from exc

    from .grids import build_optuna_distributions

    # ---- Parallelism budget (mirrors _fit_grid_search) --------------------
    if model_name == "mlp":
        _mlp_has_gpu = False
        try:
            import torch as _torch

            _mlp_has_gpu = _torch.cuda.is_available()
            try:
                _mlp_has_gpu = _mlp_has_gpu or _torch.accelerator.is_available()
            except AttributeError:
                vlog(
                    verbose,
                    "torch.accelerator not available (older torch); using cuda/mps only",
                    level="debug",
                )
            if hasattr(_torch.backends, "mps"):
                _mlp_has_gpu = _mlp_has_gpu or _torch.backends.mps.is_available()
        except ImportError:
            vlog(verbose, "torch not installed; MLP will run on CPU", level="debug")
        search_n_jobs = 1 if _mlp_has_gpu else max_cores
        model_n_jobs = 1
    elif (
        (model_name in {"xgb", "rf"} and bool(xgb_use_gpu))
        or (model_name == "svm" and bool(cuml_use_gpu))
        or (model_name == "lgbm" and bool(lgbm_use_gpu))
    ):
        search_n_jobs = 1
        model_n_jobs = 1
    elif model_name in {"rf", "xgb", "lgbm"} and max_cores > 1:
        model_n_jobs = max(1, int(max_cores**0.5))
        search_n_jobs = max(1, max_cores // model_n_jobs)
    else:
        search_n_jobs = max_cores
        model_n_jobs = 1

    # ---- Imbalance weight for XGB/RF classification (lgbm computed here too
    # for backward-compatible threading, but ignored by _build_lgbm_base_model
    # -- unsafe on device_type="cuda", applied via sample_weight at fit time
    # instead; see _lgbm_classification_sample_weight) -------------------------
    scale_pos_weight: float | None = None
    if task == "classification" and model_name in {"xgb", "rf", "lgbm"}:
        pos = int(np.sum(y_train == 1))
        neg = int(np.sum(y_train == 0))
        if pos > 0 and neg > 0:
            scale_pos_weight = neg / pos

    # ---- Helper: build sampler + OptunaSearchCV ---------------------------
    def _make_search(
        xgb_gpu: bool | None,
        cuml_gpu: bool | None,
        n_jobs: int,
        lgbm_gpu: bool | None = lgbm_use_gpu,
    ) -> OptunaSearchCV:
        """Internal helper for make search."""
        base_model, distributions = build_optuna_distributions(
            model_name=model_name,
            task=task,
            n_samples=x_train.shape[0],
            n_features=x_train.shape[1],
            xgb_use_gpu=xgb_gpu,
            cuml_use_gpu=cuml_gpu,
            lgbm_use_gpu=lgbm_gpu,
            verbose=verbose,
            model_n_jobs=model_n_jobs,
            scale_pos_weight=scale_pos_weight,
        )

        # For XGB with early stopping, n_estimators is determined by ES — not by
        # the search.  Fix it on the base model and remove from distributions so
        # Optuna only tunes the remaining hyperparameters.
        if model_name == "xgb":
            base_model.set_params(n_estimators=_XGB_ES_N_MAX)
            distributions.pop("model__n_estimators", None)
            vlog(
                verbose,
                f"XGB Optuna search: fixed n_estimators={_XGB_ES_N_MAX} (determined by ES); "
                f"removed from distributions",
                level="debug",
            )

        # Same treatment for LightGBM.
        if model_name == "lgbm":
            base_model.set_params(n_estimators=_LGBM_ES_N_MAX)
            distributions.pop("model__n_estimators", None)
            vlog(
                verbose,
                f"LGBM Optuna search: fixed n_estimators={_LGBM_ES_N_MAX} (determined by ES); "
                f"removed from distributions",
                level="debug",
            )

        pipe_copy = clone(pipe)
        pipe_copy.set_params(model=base_model)

        # Wrap pipeline with _ESPipeline for models that need early stopping
        if model_name in {"xgb", "mlp", "lgbm"}:
            pipe_copy = _ESPipeline(pipe_copy.steps, model_name=model_name)

        # Activate _ContextAwareCV for models that need inner-fold validation data
        cv_splitter = inner_cv
        if model_name in {"xgb", "mlp", "lgbm"}:
            cv_splitter = _ContextAwareCV(inner_cv, model_name)

        warnings.filterwarnings(
            "ignore", category=optuna.exceptions.ExperimentalWarning
        )
        if optuna_sampler == "cmaes":
            sampler = optuna.samplers.CmaEsSampler(
                seed=random_seed,
                n_startup_trials=optuna_n_startup_trials,
                restart_strategy="ipop",
            )
        else:
            sampler = optuna.samplers.TPESampler(
                seed=random_seed,
                n_startup_trials=optuna_n_startup_trials,
                multivariate=optuna_multivariate,
                constant_liar=True,
            )

        # Match Optuna log verbosity to the pipeline's --verbose/--debug flags.
        optuna.logging.set_verbosity(
            optuna.logging.INFO if verbose else optuna.logging.WARNING
        )
        study = optuna.create_study(direction="maximize", sampler=sampler)

        callbacks = []
        if optuna_wandb_callback:
            try:
                import wandb as _wandb
                from optuna_integration.wandb import WeightsAndBiasesCallback

                if _wandb.run is not None:
                    callbacks.append(
                        WeightsAndBiasesCallback(
                            metric_name=f"{model_name}/{scorer_name(task)}",
                            as_multirun=False,
                        )
                    )
            except Exception as exc:
                vlog(verbose, f"wandb callback setup failed: {exc}", level="debug")

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=optuna.exceptions.ExperimentalWarning
            )
            return OptunaSearchCV(
                estimator=pipe_copy,
                param_distributions=distributions,
                scoring=scorer_name(task),
                cv=cv_splitter,
                n_trials=param_grid_size,
                n_jobs=n_jobs,
                refit=True,
                error_score=np.nan,
                study=study,
                callbacks=callbacks if callbacks else None,
                verbose=0,
            )

    search = _make_search(xgb_use_gpu, cuml_use_gpu, search_n_jobs, lgbm_use_gpu)

    def _fit_search(search_obj: OptunaSearchCV) -> None:
        """Internal helper for fit search."""
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=optuna.exceptions.ExperimentalWarning
            )
            if search_n_jobs > 1:
                with parallel_backend("threading", n_jobs=search_n_jobs):
                    search_obj.fit(x_train, y_train)
                return
            search_obj.fit(x_train, y_train)

    def _fallback_to_grid(reason: str) -> tuple[Pipeline, dict[str, Any]]:
        """Fallback to the deterministic grid-search path when Optuna produces no usable trials."""
        vlog(
            verbose,
            f"Optuna search fallback to grid for model={model_name}, task={task} ({reason})",
            level="info",
        )
        fallback_xgb_use_gpu = False if model_name in {"xgb", "rf"} else xgb_use_gpu
        fallback_cuml_use_gpu = False if model_name == "svm" else cuml_use_gpu
        fallback_lgbm_use_gpu = False if model_name == "lgbm" else lgbm_use_gpu
        return _fit_grid_search(
            task=task,
            model_name=model_name,
            x_train=x_train,
            y_train=y_train,
            pipe=pipe,
            inner_cv=inner_cv,
            max_cores=max_cores,
            param_grid_size=param_grid_size,
            search_strategy="grid",
            lhs_scale_mode="auto",
            xgb_use_gpu=fallback_xgb_use_gpu,
            cuml_use_gpu=fallback_cuml_use_gpu,
            lgbm_use_gpu=fallback_lgbm_use_gpu,
            verbose=verbose,
            debug_grid_progress=False,
        )

    optuna_failure_reason: str | None = None
    _log_early_stopping_fold_sizes(model_name, inner_cv, verbose=verbose)
    try:
        _fit_search(search)
    except Exception as exc:
        # On GPU failure rebuild with a fresh study (NaN GPU trials corrupt the
        # surrogate model) and retry on CPU once.
        if (
            model_name in {"xgb", "rf"}
            and bool(xgb_use_gpu)
            and _looks_like_gpu_failure(exc)
        ):
            vlog(
                verbose,
                f"{model_name.upper()} GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            search = _make_search(False, cuml_use_gpu, max_cores)
            _fit_search(search)
        elif (
            model_name == "svm" and bool(cuml_use_gpu) and _looks_like_gpu_failure(exc)
        ):
            vlog(
                verbose,
                f"cuML SVM GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            search = _make_search(xgb_use_gpu, False, max_cores)
            _fit_search(search)
        elif (
            model_name == "lgbm" and bool(lgbm_use_gpu) and _looks_like_gpu_failure(exc)
        ):
            vlog(
                verbose,
                f"LGBM GPU training failed; retrying on CPU ({exc})",
                level="info",
            )
            search = _make_search(xgb_use_gpu, cuml_use_gpu, max_cores, False)
            _fit_search(search)
        else:
            message = str(exc)
            if (
                "No trials are completed yet" in message
                or "The value nan is not acceptable" in message
                or "failed with value np.float64(nan)" in message
            ):
                optuna_failure_reason = message
            else:
                raise

    if optuna_failure_reason is None:
        try:
            best_score = float(search.best_score_)
            if not np.isfinite(best_score):
                optuna_failure_reason = f"best_score={best_score!r}"
        except Exception as exc:
            optuna_failure_reason = str(exc)

    if optuna_failure_reason is not None:
        return _fallback_to_grid(optuna_failure_reason)

    # XGB/LightGBM/MLP: post-search refit with early stopping to find optimal n_estimators/epochs.
    tree_es_applied = False
    tree_es_n_estimators: int | None = None
    calibration_holdout: tuple[pd.DataFrame, np.ndarray] | None = None
    if model_name in {"xgb", "mlp", "lgbm"}:
        pipe_es, tree_es_applied, tree_es_n_estimators, calibration_holdout = (
            _fit_tree_early_stopping(
                model_name=model_name,
                task=task,
                fitted_pipe=search.best_estimator_,
                x_train=x_train,
                y_train=y_train,
                inner_cv=inner_cv,
                verbose=verbose,
                reserve_calibration_holdout=calibrate and model_name in {"xgb", "lgbm"},
            )
        )
        search.best_estimator_ = pipe_es

    vlog(
        verbose,
        f"Optuna search complete model={model_name}, task={task}, "
        f"sampler={optuna_sampler}, n_trials={param_grid_size}, "
        f"best_score={float(search.best_score_):.6f}",
    )

    _set_xgb_cpu_predictor_for_inference(search.best_estimator_)
    _set_lgbm_cpu_predictor_for_inference(search.best_estimator_)

    tuning_info: dict[str, Any] = {
        "best_params": sanitize_best_params(search.best_params_),
        "best_score": float(search.best_score_),
        "optuna_n_trials": int(param_grid_size),
        "optuna_sampler": optuna_sampler,
    }
    if model_name in {"xgb", "mlp", "rf", "lgbm"}:
        tuning_info["tree_es_applied"] = tree_es_applied
        tuning_info["tree_es_n_estimators"] = tree_es_n_estimators
    if calibration_holdout is not None:
        tuning_info["calibration_holdout"] = calibration_holdout

    try:
        cv_res = search.cv_results_
        vlog(verbose, f"optuna cv_results keys: {sorted(cv_res.keys())}", level="debug")
        split_keys = sorted(
            k for k in cv_res if k.startswith("split") and k.endswith("_test_score")
        )
        params = cv_res.get("params", [])
        n = len(params)
        tuning_info["cv_results"] = {
            "n_candidates": n,
            "mean_test_score": [
                float(v) if np.isfinite(float(v)) else None
                for v in cv_res.get("mean_test_score", [])
            ],
            "std_test_score": [
                float(v) if np.isfinite(float(v)) else None
                for v in cv_res.get("std_test_score", [])
            ],
            "rank_test_score": [int(v) for v in cv_res.get("rank_test_score", [])],
            "params_json": [
                json.dumps(sanitize_best_params(p), sort_keys=True, default=str)
                for p in params
            ],
            "split_keys": split_keys,
            "split_scores": {
                k: [float(v) if np.isfinite(float(v)) else None for v in cv_res[k]]
                for k in split_keys
            },
        }
    except Exception as exc:
        vlog(verbose, f"optuna cv_results extraction failed: {exc}", level="debug")

    return search.best_estimator_, tuning_info


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def fit_best_estimator(
    task: str,
    model_name: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    param_grid_size: int,
    search_strategy: str,
    lhs_scale_mode: str,
    xgb_use_gpu: bool | None,
    cuml_use_gpu: bool | None = None,
    verbose: bool = False,
    debug_grid_progress: bool = False,
    optuna_backend: bool = False,
    optuna_sampler: str = "tpe",
    optuna_n_startup_trials: int = 5,
    optuna_multivariate: bool = True,
    optuna_wandb_callback: bool = False,
    calibrate: bool = False,
    lgbm_use_gpu: bool | None = None,
    random_seed: int = RNG_SEED,
) -> tuple[Pipeline, dict[str, Any]]:
    """Tune hyperparameters via grouped inner CV and return the best estimator.

    Dispatches to one of three strategies based on model_name:
    - "elasticnet": ElasticNetCV / LogisticRegressionCV self-select alpha (C) and l1_ratio.
    - "beta": no hyperparameters; cross_validate + refit.
    - All others: GridSearchCV with grid or LHS candidates.

    calibrate (XGB and LightGBM only, via GridSearchCV/Optuna paths): reserves
    the ES-val slice as a dedicated leakage-free calibration holdout — see
    ``_fit_tree_early_stopping``'s ``reserve_calibration_holdout``.

    Returns
    -------
    (best_estimator, tuning_info)
        best_estimator is a fitted sklearn Pipeline ready for outer-fold eval.
        tuning_info contains best_params, best_score, and model-specific extras.
    """
    # ElasticNetCV / LogisticRegressionCV: self-select regularisation via built-in CV;
    # bypass GridSearchCV and pass grouped inner splits + n_jobs directly.
    if task == "regression" and model_name == "elasticnet":
        return _fit_elasticnet_cv(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            inner_cv=inner_cv,
            max_cores=max_cores,
            verbose=verbose,
            cuml_use_gpu=cuml_use_gpu,
            debug_grid_progress=debug_grid_progress,
        )

    if task == "classification" and model_name == "elasticnet":
        return _fit_logistic_elasticnet_cv(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            inner_cv=inner_cv,
            max_cores=max_cores,
            verbose=verbose,
            cuml_use_gpu=cuml_use_gpu,
            debug_grid_progress=debug_grid_progress,
        )

    # BetaRegressor: no external hyperparameters; use cross_validate for scoring.
    if task == "regression" and model_name == "beta":
        return _fit_beta(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            inner_cv=inner_cv,
            max_cores=max_cores,
            verbose=verbose,
        )

    # TabICLv2: pre-trained in-context learning model, no hyperparameters, no inner CV.
    if model_name == "tabicl":
        return _fit_tabicl(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            verbose=verbose,
        )

    # All other models: GridSearchCV or Optuna search.
    from sklearn.linear_model import LinearRegression

    pipe = Pipeline(steps=[("prep", preprocessor), ("model", LinearRegression())])

    # Optuna path: Bayesian optimisation via TPE or CMA-ES.
    # "linear" is excluded (no hyperparameters to tune).
    if optuna_backend and model_name not in {"linear"}:
        return _fit_optuna_search(
            task=task,
            model_name=model_name,
            x_train=x_train,
            y_train=y_train,
            pipe=pipe,
            inner_cv=inner_cv,
            max_cores=max_cores,
            param_grid_size=param_grid_size,
            optuna_sampler=optuna_sampler,
            optuna_n_startup_trials=optuna_n_startup_trials,
            optuna_multivariate=optuna_multivariate,
            xgb_use_gpu=xgb_use_gpu,
            cuml_use_gpu=cuml_use_gpu,
            verbose=verbose,
            optuna_wandb_callback=optuna_wandb_callback,
            calibrate=calibrate,
            lgbm_use_gpu=lgbm_use_gpu,
            random_seed=random_seed,
        )

    return _fit_grid_search(
        task=task,
        model_name=model_name,
        x_train=x_train,
        y_train=y_train,
        pipe=pipe,
        inner_cv=inner_cv,
        max_cores=max_cores,
        param_grid_size=param_grid_size,
        search_strategy=search_strategy,
        lhs_scale_mode=lhs_scale_mode,
        xgb_use_gpu=xgb_use_gpu,
        cuml_use_gpu=cuml_use_gpu,
        verbose=verbose,
        debug_grid_progress=debug_grid_progress,
        calibrate=calibrate,
        lgbm_use_gpu=lgbm_use_gpu,
        random_seed=random_seed,
    )
