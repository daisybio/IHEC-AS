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

import math
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
from ..metrics import regression_metrics
from ..utils import progress_iter, safe_json, sanitize_best_params, vlog
from .beta import BetaRegressor, _inverse_logit
from .grids import (
    build_param_candidates,
    _use_lhs_strategy,
    _ALPHA_GRID,
    _C_GRID,
    _L1_RATIO_GRID,
)
from .xgb_utils import _set_xgb_cpu_predictor_for_inference, xgb_gpu_available

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
    return (
        "balanced_accuracy"
        if task == "classification"
        else "neg_root_mean_squared_error"
    )


def metric_key(task: str) -> str:
    """Return the primary metric name for a task."""
    return "balanced_accuracy" if task == "classification" else "rmse"


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
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit ElasticNetCV — alpha and l1_ratio selection via built-in CV.

    The same grouped inner-fold splits used by all other models are passed
    directly to ElasticNetCV so fold assignments are consistent across models.
    Always sweeps _ALPHA_GRID alphas and _L1_RATIO_GRID l1_ratios.

    When ``cuml_use_gpu`` is True (or auto-detected), uses cuML ElasticNet via
    GridSearchCV over a coarser alpha grid.  cuML lacks a regularization-path
    CV variant, so individual fits are done instead.
    """
    from .cuml_utils import cuml_gpu_available

    use_gpu = bool(cuml_use_gpu) if cuml_use_gpu is not None else cuml_gpu_available()

    if use_gpu:
        return _fit_cuml_elasticnet_cv(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            inner_cv=inner_cv,
            verbose=verbose,
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
                    alphas=_ALPHA_GRID,
                    l1_ratio=_L1_RATIO_GRID,
                    n_jobs=max_cores,
                ),
            ),
        ]
    )
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
    return estimator, {
        "best_params": sanitize_best_params(
            {
                "model": model,
                "model__alpha": alpha,
                "model__l1_ratio": l1_ratio,
                "model__cv_folds": int(len(inner_cv)),
                "model__alpha_count": len(_ALPHA_GRID),
            }
        ),
        "best_score": best_score,
    }


# Coarser alpha grid for cuML ElasticNet: individual fits instead of path.
# 15 log-spaced values cover the same range as _ALPHA_GRID without the 100-point expense.
_CUML_ALPHA_GRID: np.ndarray = np.logspace(-6, 2, num=15)


def _fit_cuml_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit cuML ElasticNet via GridSearchCV over a coarser alpha/l1_ratio grid.

    cuML ElasticNet has the same alpha and l1_ratio parameters as sklearn but
    no CV-path variant.  GridSearchCV drives the search (n_jobs=1 for GPU).
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
        for a in _CUML_ALPHA_GRID
        for r in _L1_RATIO_GRID
    ]
    search = GridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=inner_cv,
        n_jobs=1,
        refit=True,
        error_score="raise",
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
                            alphas=_ALPHA_GRID,
                            l1_ratio=_L1_RATIO_GRID,
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
                f"candidates={len(_ALPHA_GRID)}, best_score={best_score:.6f}",
            )
            return cpu_est, {
                "best_params": sanitize_best_params(
                    {
                        "model": cpu_model,
                        "model__alpha": alpha,
                        "model__l1_ratio": l1_ratio,
                        "model__cv_folds": int(len(inner_cv)),
                        "model__alpha_count": len(_ALPHA_GRID),
                    }
                ),
                "best_score": best_score,
            }
        raise

    best_params = search.best_params_
    total_candidates = len(_CUML_ALPHA_GRID) * len(_L1_RATIO_GRID)
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
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit LogisticRegressionCV — C and l1_ratio selection via built-in CV.

    Mirrors the regression ElasticNetCV path: the same grouped inner-fold splits
    are passed directly so fold assignments are consistent across models, and
    n_jobs is set to max_cores.  Sweeps _C_GRID C values and _L1_RATIO_GRID
    l1_ratios (C = 1 / alpha, reciprocal of _ALPHA_GRID).

    When ``cuml_use_gpu`` is True (or auto-detected), uses cuML LogisticRegression
    via GridSearchCV over a coarser C/l1_ratio grid.
    """
    from .cuml_utils import cuml_gpu_available

    use_gpu = bool(cuml_use_gpu) if cuml_use_gpu is not None else cuml_gpu_available()

    if use_gpu:
        return _fit_cuml_logistic_elasticnet_cv(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            inner_cv=inner_cv,
            verbose=verbose,
        )

    from sklearn.linear_model import LogisticRegressionCV

    estimator = Pipeline(
        steps=[
            ("prep", preprocessor),
            (
                "model",
                LogisticRegressionCV(
                    Cs=_C_GRID,
                    l1_ratios=_L1_RATIO_GRID,
                    solver="saga",
                    max_iter=20000,
                    class_weight="balanced",
                    random_state=RNG_SEED,
                    cv=inner_cv,
                    n_jobs=max_cores,
                    use_legacy_attributes=False,
                ),
            ),
        ]
    )
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
        except Exception:
            pass

    vlog(
        verbose,
        f"Search complete model=elasticnet, task={task}, strategy=logisticelasticnetcv, "
        f"candidates={n_cs if n_cs > 0 else len(_C_GRID)}, best_score={best_score:.6f}",
    )
    return estimator, {
        "best_params": sanitize_best_params(
            {
                "model": model,
                "model__C": best_c,
                "model__l1_ratio": l1_ratio,
                "model__cv_folds": int(len(inner_cv)),
                "model__Cs_count": len(_C_GRID),
            }
        ),
        "best_score": best_score,
    }


# Coarser C grid for cuML LogisticRegression: individual fits instead of path.
# 15 log-spaced values on the same range as _C_GRID (= 1 / _ALPHA_GRID).
_CUML_C_GRID: np.ndarray = np.logspace(-6, 2, num=15)[::-1]  # high→low (descending C)


def _fit_cuml_logistic_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit cuML LogisticRegression via GridSearchCV over a coarser C/l1_ratio grid.

    cuML LogisticRegression supports penalty="elasticnet", class_weight="balanced",
    C, and l1_ratio with the same semantics as sklearn.  No CV-path variant exists,
    so GridSearchCV drives the search (n_jobs=1 for GPU).
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
        for c in _CUML_C_GRID
        for r in _L1_RATIO_GRID
    ]
    search = GridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=inner_cv,
        n_jobs=1,
        refit=True,
        error_score="raise",
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
            from sklearn.linear_model import LogisticRegressionCV as _LogisticRegressionCV

            cpu_est = Pipeline(
                steps=[
                    ("prep", preprocessor),
                    (
                        "model",
                        _LogisticRegressionCV(
                            Cs=_C_GRID,
                            l1_ratios=_L1_RATIO_GRID,
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
            l1_ratio = float(np.asarray(getattr(cpu_model, "l1_ratio_", np.nan)).ravel()[0])
            best_score = float("nan")
            raw_scores = getattr(cpu_model, "scores_", None)
            if raw_scores is not None:
                try:
                    scores_arr = (
                        next(iter(raw_scores.values()))
                        if isinstance(raw_scores, dict)
                        else raw_scores
                    )
                    best_score = float(np.max(np.mean(np.asarray(scores_arr, dtype=float), axis=0)))
                except Exception:
                    pass
            vlog(
                verbose,
                f"Search complete model=elasticnet, task={task}, strategy=logisticelasticnetcv_cpu_fallback, "
                f"candidates={len(_C_GRID)}, best_score={best_score:.6f}",
            )
            return cpu_est, {
                "best_params": sanitize_best_params(
                    {
                        "model": cpu_model,
                        "model__C": best_c,
                        "model__l1_ratio": l1_ratio,
                        "model__cv_folds": int(len(inner_cv)),
                        "model__Cs_count": len(_C_GRID),
                    }
                ),
                "best_score": best_score,
            }
        raise

    best_params = search.best_params_
    total_candidates = len(_CUML_C_GRID) * len(_L1_RATIO_GRID)
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


def _looks_like_gpu_failure(exc: Exception) -> bool:
    """Return True when an exception message suggests a CUDA/GPU problem."""
    msg = str(exc).lower()
    return any(
        token in msg
        for token in ("cuda", "gpu", "out of memory", "oom", "device", "cublas", "nccl")
    )


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
) -> GridSearchCV:
    """Construct a GridSearchCV object with the appropriate candidate list."""
    scale_pos_weight = None
    if model_name == "rf" and task == "classification":
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
        verbose=verbose,
        model_n_jobs=model_n_jobs,
        scale_pos_weight=scale_pos_weight,
    )
    return _ProgressGridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=inner_cv,
        n_jobs=search_n_jobs,
        pre_dispatch=search_n_jobs,
        refit=True,
        error_score="raise",
        progress_enabled=debug_grid_progress,
        progress_desc=f"grid {model_name} [{task}]",
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
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit all non-ElasticNet, non-Beta models via GridSearchCV.

    Handles:
    - Parallelism: limits XGBoost/cuML GPU jobs to 1 to avoid VRAM contention.
    - GPU retry: on CUDA failure, automatically falls back to CPU.
    - XGBoost post-search refit: after tuning with a reduced n_estimators
      budget, the final estimator is re-trained with a larger tree budget.
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
                pass
            if hasattr(_torch.backends, "mps"):
                _mlp_has_gpu = _mlp_has_gpu or _torch.backends.mps.is_available()
        except ImportError:
            pass
        # Single search job when a GPU is available to prevent VRAM contention
        # across parallel GridSearchCV workers.
        search_n_jobs = 1 if _mlp_has_gpu else max_cores
        model_n_jobs = 1
    elif (model_name in {"xgb", "rf"} and bool(xgb_use_gpu)) or (
        model_name == "svm" and bool(cuml_use_gpu)
    ):
        search_n_jobs = 1
        model_n_jobs = 1
    elif model_name in {"rf", "xgb"} and max_cores > 1:
        # Split budget between search-level and model-level parallelism.
        # RF uses joblib tree-building threads; XGB CPU uses its own thread pool.
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
        verbose=verbose,
        model_n_jobs=model_n_jobs,
        debug_grid_progress=debug_grid_progress,
    )

    def _fit_search(search_obj: GridSearchCV) -> None:
        # sklearn 1.8 internally converts C=np.inf → penalty=None and then
        # warns about it, even though C=np.inf is their own recommended API.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Setting penalty=None will ignore the C and l1_ratio parameters",
                category=UserWarning,
            )
            if search_n_jobs > 1:
                # Thread backend avoids occasional loky worker-stop warnings.
                with parallel_backend("threading", n_jobs=search_n_jobs):
                    search_obj.fit(x_train, y_train)
                return
            search_obj.fit(x_train, y_train)

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
                verbose=verbose,
                debug_grid_progress=debug_grid_progress,
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
                verbose=verbose,
                debug_grid_progress=debug_grid_progress,
            )
            _fit_search(search)
        else:
            raise

    # XGBoost: tune with reduced n_estimators for speed, then refit with full budget.
    xgb_refit_applied = False
    xgb_final_n_estimators: int | None = None
    if model_name == "xgb":
        try:
            best_params = search.best_params_.copy()
            n_est_full = 300 if x_train.shape[0] >= 5000 else 180
            best_params["model__n_estimators"] = n_est_full
            if bool(xgb_use_gpu):
                # Refit on CPU to prevent post-search GPU OOM.
                best_params["model__device"] = "cpu"
            pipe_final = clone(pipe)
            pipe_final.set_params(**best_params)
            pipe_final.fit(x_train, y_train)
            search.best_estimator_ = pipe_final
            xgb_refit_applied = True
            xgb_final_n_estimators = n_est_full
        except Exception as exc:
            # Keep the tuned estimator if the full-budget refit fails.
            vlog(
                verbose,
                f"XGB full-refit skipped; keeping tuned estimator ({exc})",
                level="info",
            )

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

    tuning_info: dict[str, Any] = {
        "best_params": sanitize_best_params(search.best_params_),
        "best_score": float(search.best_score_),
    }
    if model_name == "xgb":
        tuning_info["xgb_refit_applied"] = xgb_refit_applied
        tuning_info["xgb_final_n_estimators"] = xgb_final_n_estimators

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
) -> tuple[Pipeline, dict[str, Any]]:
    """Tune hyperparameters via grouped inner CV and return the best estimator.

    Dispatches to one of three strategies based on model_name:
    - "elasticnet": ElasticNetCV / LogisticRegressionCV self-select alpha (C) and l1_ratio.
    - "beta": no hyperparameters; cross_validate + refit.
    - All others: GridSearchCV with grid or LHS candidates.

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

    # All other models: GridSearchCV with candidate grid or LHS samples.
    from sklearn.linear_model import LinearRegression

    pipe = Pipeline(steps=[("prep", preprocessor), ("model", LinearRegression())])
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
    )
