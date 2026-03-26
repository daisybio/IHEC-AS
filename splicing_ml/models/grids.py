from __future__ import annotations

"""Hyperparameter candidate builders for GridSearchCV.

Public entry point: ``build_param_candidates`` — a single function that
dispatches to grid or Latin Hypercube sampling based on the requested
search strategy and model type.

The previous three-function API (``choose_param_grid``,
``choose_param_distributions``, ``choose_param_lhs_candidates``) is reduced to
two internal implementations because ``choose_param_distributions`` was dead
code (the pipeline always uses GridSearchCV with list-of-dicts, never
RandomizedSearchCV with distributions).
"""

import math
from typing import Any

import numpy as np

from ..config import RNG_SEED
from ..utils import vlog

__all__ = [
    "build_param_candidates",
    # Legacy names kept for backward compatibility via the modeling shim.
    "choose_param_grid",
    "choose_param_lhs_candidates",
    "_lhs_unit",
    "_map_with_scale",
    "_axis_count_from_budget",
    "_pick_evenly_spaced",
]


# ---------------------------------------------------------------------------
# Low-level sampling helpers
# ---------------------------------------------------------------------------


def _axis_count_from_budget(grid_size: int, n_dims: int) -> int:
    """Return per-axis candidate count to target an overall grid-size budget."""
    budget = max(1, int(grid_size))
    dims = max(1, int(n_dims))
    return max(1, int(round(budget ** (1.0 / float(dims)))))


def _pick_evenly_spaced(values: list[Any], k: int) -> list[Any]:
    """Pick k values spread evenly across a candidate list (order preserved)."""
    if not values:
        return []
    if k >= len(values):
        return list(values)
    idx = np.linspace(0, len(values) - 1, num=max(1, k))
    picked_idx = sorted({int(round(i)) for i in idx})
    return [values[i] for i in picked_idx]


def _lhs_unit(n_points: int, n_dims: int, seed: int = RNG_SEED) -> np.ndarray:
    """Generate Latin Hypercube sample points in [0, 1]^d.

    Falls back to uniform random sampling when scipy.stats.qmc is unavailable.
    """
    try:
        from scipy.stats import qmc

        sampler = qmc.LatinHypercube(d=max(1, int(n_dims)), seed=seed)
        return sampler.random(n=max(1, int(n_points)))
    except Exception:
        rng = np.random.default_rng(seed)
        return rng.random((max(1, int(n_points)), max(1, int(n_dims))))


def _map_log10(u: float, low: float, high: float) -> float:
    """Map a unit [0, 1] value to [low, high] on log10 scale."""
    return float(10.0 ** (np.log10(low) + u * (np.log10(high) - np.log10(low))))


def _map_linear(u: float, low: float, high: float) -> float:
    """Map a unit [0, 1] value to [low, high] on linear scale."""
    return float(low + u * (high - low))


def _map_with_scale(
    u: float,
    low: float,
    high: float,
    scale_mode: str,
    default_scale: str,
) -> float:
    """Map unit value to [low, high] using the configured or default scale.

    Parameters
    ----------
    scale_mode : {"auto", "log", "linear"}
        "auto" defers to ``default_scale``; "log" and "linear" override it.
    default_scale : {"log", "linear"}
        Scale used when ``scale_mode`` is "auto".
    """
    chosen = default_scale if scale_mode == "auto" else scale_mode
    if chosen == "log" and low > 0 and high > 0:
        return _map_log10(u, low, high)
    return _map_linear(u, low, high)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def build_param_candidates(
    model_name: str,
    task: str,
    n_samples: int,
    n_features: int,
    budget: int,
    strategy: str = "grid",
    lhs_scale_mode: str = "auto",
    xgb_use_gpu: bool | None = None,
    cv: list[tuple[Any, Any]] | None = None,
    verbose: bool = False,
    model_n_jobs: int = 1,
) -> list[dict[str, Any]]:
    """Return a GridSearchCV-compatible candidate list for the given model.

    Dispatches to LHS-based sampling (``choose_param_lhs_candidates``) when
    the effective strategy calls for it, or to the fixed discrete grid
    (``choose_param_grid``) otherwise.

    Parameters
    ----------
    model_name : str
        One of: "linear", "lasso", "beta", "rf", "xgb".
    task : str
        "regression" or "classification".
    n_samples, n_features : int
        Dataset dimensions used to scale grid sizes heuristically.
    budget : int
        Target number of hyperparameter candidates.
    strategy : {"grid", "random", "hybrid"}
        "random" and "hybrid" (for rf/xgb) use LHS; "grid" uses fixed pools.
    lhs_scale_mode : {"auto", "log", "linear"}
        Scale mapping for LHS candidates (passed through to LHS builder).
    xgb_use_gpu : bool or None
        GPU flag forwarded to XGBoost candidate builders.
    cv : list of (train_idx, val_idx) tuples, optional
        Passed through for informational logging only; not consumed here.
    verbose : bool
        Enable debug logging.
    """
    use_lhs = _use_lhs_strategy(strategy, model_name)
    if use_lhs:
        return choose_param_lhs_candidates(
            model_name=model_name,
            task=task,
            n_samples=n_samples,
            n_features=n_features,
            budget=budget,
            xgb_use_gpu=xgb_use_gpu,
            scale_mode=lhs_scale_mode,
            cv=cv,
            verbose=verbose,
            model_n_jobs=model_n_jobs,
        )
    return choose_param_grid(
        model_name=model_name,
        task=task,
        n_samples=n_samples,
        n_features=n_features,
        max_cores=1,  # max_cores is only used for logging in choose_param_grid
        grid_size=budget,
        xgb_use_gpu=xgb_use_gpu,
        cv=cv,
        verbose=verbose,
        model_n_jobs=model_n_jobs,
    )


def _use_lhs_strategy(search_strategy: str, model_name: str) -> bool:
    """Return True when LHS sampling should be used instead of a fixed grid."""
    if search_strategy == "random":
        return True
    if search_strategy == "grid":
        return False
    # Hybrid: use LHS for high-dimensional tree model spaces.
    return model_name in {"rf", "xgb"}


# ---------------------------------------------------------------------------
# Fixed-grid implementation
# ---------------------------------------------------------------------------


def choose_param_grid(
    model_name: str,
    task: str,
    n_samples: int,
    n_features: int,
    max_cores: int,
    grid_size: int,
    xgb_use_gpu: bool | None = None,
    cv: list[tuple[Any, Any]] | None = None,
    verbose: bool = False,
    model_n_jobs: int = 1,
) -> list[dict[str, Any]]:
    """Build heuristic fixed-grid parameter candidates for GridSearchCV.

    Each returned dict maps pipeline parameter names to lists of candidate
    values. The "model" key holds the estimator instance(s) to try.
    """
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.linear_model import (
        Lasso,
        LassoCV,
        LinearRegression,
        LogisticRegression,
        LogisticRegressionCV,
        Ridge,
    )

    from .beta import BetaRegressor
    from .xgb_utils import xgb_gpu_available

    small = n_samples < 5000
    high_dim = n_features > n_samples
    vlog(
        verbose,
        f"Choosing param grid for model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}, "
        f"max_cores={max_cores}, grid_size={grid_size}",
    )

    if model_name == "linear" and task == "classification":
        c_base = max(0.01, min(10.0, n_samples / 1000.0))
        c_candidates = [
            c_base / 30.0,
            c_base / 10.0,
            c_base,
            c_base * 10.0,
            c_base * 30.0,
        ]
        c_grid = _pick_evenly_spaced(c_candidates, max(1, min(grid_size, 5)))
        return [
            {
                "model": [
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=RNG_SEED,
                    )
                ],
                "model__C": c_grid,
                "model__solver": ["lbfgs"],
            }
        ]

    if model_name == "lasso" and task == "classification":
        # Lasso-equivalent for classification: L1-regularised logistic regression.
        # LogisticRegressionCV sweeps 100 C values internally (one fit), analogous
        # to LassoCV for regression.
        return [
            {
                "model": [
                    LogisticRegressionCV(
                        Cs=100,
                        l1_ratios=(1,),
                        solver="saga",
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=RNG_SEED,
                        cv=5,
                        n_jobs=1,
                    )
                ],
            }
        ]

    if model_name == "linear" and task == "regression":
        if high_dim:
            ridge_alpha = _pick_evenly_spaced(
                [0.01, 0.1, 1.0, 10.0, 100.0], max(1, min(grid_size, 5))
            )
            lasso_alpha = _pick_evenly_spaced(
                [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1],
                max(1, min(grid_size, 8)),
            )
            return [
                {
                    "model": [Ridge(random_state=RNG_SEED)],
                    "model__alpha": ridge_alpha,
                },
                {
                    "model": [Lasso(random_state=RNG_SEED, max_iter=10000)],
                    "model__alpha": lasso_alpha,
                },
            ]
        return [{"model": [LinearRegression(n_jobs=1)]}]

    if model_name == "lasso" and task == "regression":
        return [
            {
                "model": [
                    LassoCV(
                        random_state=RNG_SEED,
                        max_iter=20000,
                        cv=5,
                        alphas=np.logspace(-6, -2, num=100),
                    )
                ],
            }
        ]

    if model_name == "beta":
        if task != "regression":
            raise ValueError("beta model supports regression task only")
        return [{"model": [BetaRegressor()]}]

    if model_name == "rf":
        depth_base = max(3, int(np.log2(max(n_samples, 2))))
        axis_count = _axis_count_from_budget(grid_size, n_dims=4)
        n_estimators_pool = [80, 120, 180, 250] if small else [120, 200, 320, 480]
        depth_pool = [
            max(2, depth_base - 2),
            max(2, depth_base - 1),
            depth_base,
            depth_base + 1,
            depth_base + 2,
        ]
        min_leaf_pool = [1, 2, 4, 8]
        max_features_pool = ["sqrt", "log2"]
        n_estimators = _pick_evenly_spaced(n_estimators_pool, axis_count)
        max_depth = _pick_evenly_spaced(depth_pool, axis_count)
        min_leaf = _pick_evenly_spaced(min_leaf_pool, axis_count)
        max_features = _pick_evenly_spaced(max_features_pool, axis_count)
        if task == "classification":
            return [
                {
                    "model": [
                        RandomForestClassifier(
                            random_state=RNG_SEED,
                            n_jobs=model_n_jobs,
                            class_weight="balanced",
                        )
                    ],
                    "model__n_estimators": n_estimators,
                    "model__max_depth": max_depth,
                    "model__max_features": max_features,
                    "model__min_samples_leaf": min_leaf,
                }
            ]
        return [
            {
                "model": [RandomForestRegressor(random_state=RNG_SEED, n_jobs=model_n_jobs)],
                "model__n_estimators": n_estimators,
                "model__max_depth": max_depth,
                "model__max_features": max_features,
                "model__min_samples_leaf": min_leaf,
            }
        ]

    if model_name == "xgb":
        try:
            from xgboost import XGBClassifier, XGBRegressor
        except Exception as exc:
            raise RuntimeError("xgboost requested but not installed") from exc

        use_gpu = bool(xgb_use_gpu) if xgb_use_gpu is not None else xgb_gpu_available()
        common = {
            "random_state": RNG_SEED,
            "n_jobs": 1,
            "tree_method": "hist",
            "device": "cuda" if use_gpu else "cpu",
        }
        n_points = max(1, int(grid_size))
        depth_low, depth_high = (2, 8) if small else (3, 10)
        n_est_low, n_est_high = (60, 200) if small else (100, 500)
        unit = _lhs_unit(n_points=n_points, n_dims=5, seed=RNG_SEED)
        base_model = (
            XGBClassifier(eval_metric="logloss", **common)
            if task == "classification"
            else XGBRegressor(eval_metric="rmse", **common)
        )
        out: list[dict[str, list[Any]]] = []
        for row in unit:
            out.append(
                {
                    "model": [base_model],
                    "model__n_estimators": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[0]),
                                    n_est_low,
                                    n_est_high,
                                    "linear",
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__max_depth": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[1]),
                                    depth_low,
                                    depth_high,
                                    "linear",
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__learning_rate": [
                        _map_with_scale(
                            float(row[2]), 0.01, 0.3, "log", default_scale="log"
                        )
                    ],
                    "model__subsample": [
                        _map_with_scale(
                            float(row[3]), 0.6, 1.0, "linear", default_scale="linear"
                        )
                    ],
                    "model__colsample_bytree": [
                        _map_with_scale(
                            float(row[4]), 0.6, 1.0, "linear", default_scale="linear"
                        )
                    ],
                }
            )
        return out

    raise ValueError(f"Unsupported model name: {model_name}")


# ---------------------------------------------------------------------------
# LHS-based implementation
# ---------------------------------------------------------------------------


def choose_param_lhs_candidates(
    model_name: str,
    task: str,
    n_samples: int,
    n_features: int,
    budget: int,
    xgb_use_gpu: bool | None = None,
    scale_mode: str = "auto",
    cv: list[tuple[Any, Any]] | None = None,
    verbose: bool = False,
    model_n_jobs: int = 1,
) -> list[dict[str, list[Any]]]:
    """Build a space-filling LHS candidate set for GridSearchCV.

    Returns a list of singleton dicts compatible with
    ``GridSearchCV(param_grid=<result>)``.
    """
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.linear_model import (
        Lasso,
        LassoCV,
        LinearRegression,
        LogisticRegression,
        LogisticRegressionCV,
        Ridge,
    )

    from .beta import BetaRegressor
    from .xgb_utils import xgb_gpu_available

    n_points = max(1, int(budget))
    small = n_samples < 5000
    high_dim = n_features > n_samples
    vlog(
        verbose,
        f"Choosing LHS candidate space model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}, "
        f"budget={n_points}, scale_mode={scale_mode}",
    )

    if model_name == "linear" and task == "classification":
        c_base = max(0.01, min(10.0, n_samples / 1000.0))
        c_low, c_high = c_base / 100.0, c_base * 100.0
        unit = _lhs_unit(n_points=n_points, n_dims=1, seed=RNG_SEED)
        model = LogisticRegression(
            max_iter=5000,
            class_weight="balanced",
            random_state=RNG_SEED,
        )
        return [
            {
                "model": [model],
                "model__C": [
                    _map_with_scale(
                        float(row[0]), c_low, c_high, scale_mode, default_scale="log"
                    )
                ],
                "model__solver": ["lbfgs"],
            }
            for row in unit
        ]

    if model_name == "lasso" and task == "classification":
        # LogisticRegressionCV sweeps 100 C values internally (one fit), analogous
        # to LassoCV for regression — no LHS sampling needed.
        return [
            {
                "model": [
                    LogisticRegressionCV(
                        Cs=100,
                        l1_ratios=(1,),
                        solver="saga",
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=RNG_SEED,
                        cv=5,
                        n_jobs=1,
                    )
                ],
            }
        ]

    if model_name == "linear" and task == "regression":
        if not high_dim:
            return [{"model": [LinearRegression(n_jobs=1)]}]
        unit = _lhs_unit(n_points=n_points, n_dims=2, seed=RNG_SEED)
        out: list[dict[str, list[Any]]] = []
        for row in unit:
            if float(row[0]) < 0.5:
                out.append(
                    {
                        "model": [Ridge(random_state=RNG_SEED)],
                        "model__alpha": [
                            _map_with_scale(
                                float(row[1]),
                                1e-5,
                                1e2,
                                scale_mode,
                                default_scale="log",
                            )
                        ],
                    }
                )
            else:
                out.append(
                    {
                        "model": [Lasso(random_state=RNG_SEED, max_iter=10000)],
                        "model__alpha": [
                            _map_with_scale(
                                float(row[1]),
                                1e-6,
                                1e-1,
                                scale_mode,
                                default_scale="log",
                            )
                        ],
                    }
                )
        return out

    if model_name == "lasso" and task == "regression":
        return [
            {
                "model": [
                    LassoCV(
                        random_state=RNG_SEED,
                        max_iter=20000,
                        cv=5,
                        alphas=np.logspace(-6, -2, num=100),
                    )
                ],
            }
        ]

    if model_name == "beta":
        if task != "regression":
            raise ValueError("beta model supports regression task only")
        return [{"model": [BetaRegressor()]}]

    if model_name == "rf":
        depth_base = max(3, int(np.log2(max(n_samples, 2))))
        depth_low, depth_high = max(2, depth_base - 3), depth_base + 3
        n_est_low, n_est_high = (80, 350) if small else (120, 640)
        base_model = (
            RandomForestClassifier(
                random_state=RNG_SEED, n_jobs=model_n_jobs, class_weight="balanced"
            )
            if task == "classification"
            else RandomForestRegressor(random_state=RNG_SEED, n_jobs=model_n_jobs)
        )
        unit = _lhs_unit(n_points=n_points, n_dims=4, seed=RNG_SEED)
        out: list[dict[str, list[Any]]] = []
        for row in unit:
            out.append(
                {
                    "model": [base_model],
                    "model__n_estimators": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[0]),
                                    n_est_low,
                                    n_est_high,
                                    scale_mode,
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__max_depth": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[1]),
                                    depth_low,
                                    depth_high,
                                    scale_mode,
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__min_samples_leaf": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[2]),
                                    1.0,
                                    12.0,
                                    scale_mode,
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__max_features": ["sqrt" if float(row[3]) < 0.5 else "log2"],
                }
            )
        return out

    if model_name == "xgb":
        try:
            from xgboost import XGBClassifier, XGBRegressor
        except Exception as exc:
            raise RuntimeError("xgboost requested but not installed") from exc

        use_gpu = bool(xgb_use_gpu) if xgb_use_gpu is not None else xgb_gpu_available()
        common = {
            "random_state": RNG_SEED,
            "n_jobs": 1,
            "tree_method": "hist",
            "device": "cuda" if use_gpu else "cpu",
        }
        base_model = (
            XGBClassifier(eval_metric="logloss", **common)
            if task == "classification"
            else XGBRegressor(eval_metric="rmse", **common)
        )
        depth_low, depth_high = (2, 8) if small else (3, 10)
        n_est_low, n_est_high = (60, 200) if small else (100, 500)
        unit = _lhs_unit(n_points=n_points, n_dims=5, seed=RNG_SEED)
        out: list[dict[str, list[Any]]] = []
        for row in unit:
            out.append(
                {
                    "model": [base_model],
                    "model__n_estimators": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[0]),
                                    n_est_low,
                                    n_est_high,
                                    scale_mode,
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__max_depth": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[1]),
                                    depth_low,
                                    depth_high,
                                    scale_mode,
                                    default_scale="linear",
                                )
                            )
                        )
                    ],
                    "model__learning_rate": [
                        _map_with_scale(
                            float(row[2]), 0.01, 0.3, scale_mode, default_scale="log"
                        )
                    ],
                    "model__subsample": [
                        _map_with_scale(
                            float(row[3]), 0.6, 1.0, scale_mode, default_scale="linear"
                        )
                    ],
                    "model__colsample_bytree": [
                        _map_with_scale(
                            float(row[4]), 0.6, 1.0, scale_mode, default_scale="linear"
                        )
                    ],
                }
            )
        return out

    raise ValueError(f"Unsupported model name: {model_name}")
