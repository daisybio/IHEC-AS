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
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.pipeline import Pipeline

from ..config import RNG_SEED
from ..metrics import regression_metrics
from ..utils import safe_json, sanitize_best_params, vlog
from .beta import BetaRegressor, _inverse_logit
from .grids import build_param_candidates, _use_lhs_strategy, _ALPHA_GRID, _C_GRID, _L1_RATIO_GRID
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


def _fit_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit ElasticNetCV — alpha and l1_ratio selection via built-in CV.

    The same grouped inner-fold splits used by all other models are passed
    directly to ElasticNetCV so fold assignments are consistent across models.
    Always sweeps _ALPHA_GRID alphas and _L1_RATIO_GRID l1_ratios.
    """
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


def _fit_logistic_elasticnet_cv(
    task: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    preprocessor: ColumnTransformer,
    inner_cv: list[tuple[np.ndarray, np.ndarray]],
    max_cores: int,
    verbose: bool,
) -> tuple[Pipeline, dict[str, Any]]:
    """Fit LogisticRegressionCV — C and l1_ratio selection via built-in CV.

    Mirrors the regression ElasticNetCV path: the same grouped inner-fold splits
    are passed directly so fold assignments are consistent across models, and
    n_jobs is set to max_cores.  Sweeps _C_GRID C values and _L1_RATIO_GRID
    l1_ratios (C = 1 / alpha, reciprocal of _ALPHA_GRID).
    """
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
        cuml_use_gpu=cuml_use_gpu,
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
    cuml_use_gpu: bool | None,
    verbose: bool,
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
    elif (model_name == "xgb" and bool(xgb_use_gpu)) or (
        model_name == "svm" and bool(cuml_use_gpu)
    ):
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
        cuml_use_gpu=cuml_use_gpu,
        verbose=verbose,
        model_n_jobs=model_n_jobs,
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
                cuml_use_gpu=cuml_use_gpu,
                verbose=verbose,
            )
            _fit_search(search)
        elif model_name == "svm" and bool(cuml_use_gpu) and _looks_like_gpu_failure(exc):
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
    )
