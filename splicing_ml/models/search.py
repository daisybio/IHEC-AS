from __future__ import annotations

"""Hyperparameter tuning and best-estimator selection.

Public entry point: ``fit_best_estimator`` — dispatches to one of three
strategy helpers based on model type:

- ``_fit_lasso_cv``   : LassoCV self-selects alpha via its own internal CV.
- ``_fit_beta``       : BetaRegressor has no hyperparameters; uses
                        ``sklearn.model_selection.cross_validate`` to score
                        then refits on the full training fold.
- ``_fit_grid_search``: All other models use GridSearchCV (or LHS-based
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
from sklearn.linear_model import LassoCV
from sklearn.exceptions import FitFailedWarning
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.pipeline import Pipeline

from ..config import RNG_SEED
from ..metrics import regression_metrics
from ..utils import safe_json, sanitize_best_params, vlog
from .beta import BetaRegressor, _inverse_logit
from .grids import build_param_candidates, _use_lhs_strategy
from .xgb_utils import _set_xgb_cpu_predictor_for_inference, xgb_gpu_available

__all__ = ["fit_best_estimator"]


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


def _fit_lasso_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    param_grid_size: int,
    max_cores: int,
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit LassoCV — which handles alpha selection via its own built-in CV.

    The same grouped inner-fold splits used by all other models are passed
    directly to LassoCV so fold assignments are consistent across models.
    """
    alpha_count = max(2, int(param_grid_size))
    estimator = Pipeline(
        steps=[
            ("prep", preprocessor),
            (
                "model",
                LassoCV(
                    random_state=RNG_SEED,
                    max_iter=20000,
                    cv=inner_cv,
                    alphas=np.logspace(-6, -2, num=alpha_count),
                    n_jobs=max_cores,
                ),
            ),
        ]
    )
    estimator.fit(x_train, y_train)
    model = estimator.named_steps["model"]

    # Extract cross-validated RMSE from the MSE path produced by LassoCV.
    mse_path = np.asarray(getattr(model, "mse_path_", np.array([])))
    best_score = float("nan")
    if mse_path.size:
        best_mse = float(np.min(np.mean(mse_path, axis=1)))
        best_score = -float(np.sqrt(best_mse))

    alpha = float(getattr(model, "alpha_", np.nan))
    n_alphas = int(np.asarray(getattr(model, "alphas_", np.array([]))).size)
    vlog(
        verbose,
        f"Search complete model=lasso, task={task}, strategy=lassocv, "
        f"candidates={n_alphas if n_alphas > 0 else 1}, best_score={best_score:.6f}",
    )
    return estimator, {
        "best_params": sanitize_best_params(
            {
                "model": model,
                "model__alpha": alpha,
                "model__cv_folds": int(len(inner_cv)),
                "model__alpha_count": int(alpha_count),
            }
        ),
        "best_score": best_score,
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
    verbose: bool,
    model_n_jobs: int = 1,
) -> GridSearchCV:
    """Construct a GridSearchCV object with the appropriate candidate list."""
    candidates = build_param_candidates(
        model_name=model_name,
        task=task,
        n_samples=x_train.shape[0],
        n_features=x_train.shape[1],
        budget=param_grid_size,
        strategy=search_strategy,
        lhs_scale_mode=lhs_scale_mode,
        xgb_use_gpu=xgb_use_gpu,
        verbose=verbose,
        model_n_jobs=model_n_jobs,
    )
    return GridSearchCV(
        estimator=pipe,
        param_grid=candidates,
        scoring=scorer_name(task),
        cv=inner_cv,
        n_jobs=search_n_jobs,
        pre_dispatch=search_n_jobs,
        refit=True,
        error_score="raise",
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
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit all non-Lasso, non-Beta models via GridSearchCV.

    Handles:
    - Parallelism: limits XGBoost GPU jobs to 1 to avoid VRAM contention.
    - GPU retry: on CUDA failure, automatically falls back to CPU.
    - XGBoost post-search refit: after tuning with a reduced n_estimators
      budget, the final estimator is re-trained with a larger tree budget.
    """
    # Limit concurrent GPU jobs to avoid VRAM contention.
    # For RF, split cores between search-level and tree-level parallelism:
    # model_n_jobs gets sqrt(max_cores) cores per forest; search_n_jobs gets the rest.
    if model_name == "xgb" and bool(xgb_use_gpu):
        search_n_jobs = 1
        model_n_jobs = 1
    elif model_name in {"rf", "xgb"} and max_cores > 1:
        # Split budget between search-level and model-level parallelism.
        # RF uses joblib tree-building threads; XGB CPU uses its own thread pool.
        model_n_jobs = max(1, int(max_cores ** 0.5))
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
        verbose=verbose,
        model_n_jobs=model_n_jobs,
    )

    def _fit_search(search_obj: GridSearchCV) -> None:
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
        if model_name == "xgb" and bool(xgb_use_gpu) and _looks_like_gpu_failure(exc):
            vlog(
                verbose,
                f"XGB GPU training failed; retrying on CPU ({exc})",
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
                verbose=verbose,
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
    verbose: bool = False,
) -> tuple[Pipeline, dict[str, Any]]:
    """Tune hyperparameters via grouped inner CV and return the best estimator.

    Dispatches to one of three strategies based on model_name:
    - "lasso" (regression): LassoCV self-selects alpha.
    - "beta": no hyperparameters; cross_validate + refit.
    - All others: GridSearchCV with grid or LHS candidates.

    Returns
    -------
    (best_estimator, tuning_info)
        best_estimator is a fitted sklearn Pipeline ready for outer-fold eval.
        tuning_info contains best_params, best_score, and model-specific extras.
    """
    # LassoCV: has its own CV-based alpha selection; bypass GridSearchCV.
    if task == "regression" and model_name == "lasso":
        return _fit_lasso_cv(
            task=task,
            x_train=x_train,
            y_train=y_train,
            preprocessor=preprocessor,
            inner_cv=inner_cv,
            param_grid_size=param_grid_size,
            max_cores=max_cores,
            verbose=verbose,
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
        verbose=verbose,
    )
