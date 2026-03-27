from __future__ import annotations

"""Hyperparameter candidate builders for GridSearchCV.

Public entry point: ``build_param_candidates`` — a single function that
dispatches to grid or Latin Hypercube sampling based on the requested
search strategy and model type.

Only "linear", "rf", "xgb", "svm", and "mlp" models are handled here.
"elasticnet" and "beta" bypass GridSearchCV entirely in
``search.fit_best_estimator`` and are never dispatched to these builders.
"""

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
    # Shared regularisation grids (alpha and C are reciprocals of each other).
    "_ALPHA_GRID",
    "_C_GRID",
    "_L1_RATIO_GRID",
]

# Regularisation grids used by search._fit_elasticnet_cv and
# search._fit_logistic_elasticnet_cv.  C = 1 / alpha.
_ALPHA_GRID: np.ndarray = np.logspace(-6, 2, num=100)
_C_GRID: np.ndarray = 1.0 / _ALPHA_GRID
_L1_RATIO_GRID: list[float] = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]


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
    cuml_use_gpu: bool | None = None,
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
        One of: "linear", "rf", "xgb", "svm", "mlp".
    task : str
        "regression" or "classification".
    n_samples, n_features : int
        Dataset dimensions used to scale grid sizes heuristically.
    budget : int
        Target number of hyperparameter candidates.
    strategy : {"grid", "random", "hybrid"}
        "random" and "hybrid" (for rf/xgb/svm/mlp) use LHS; "grid" uses fixed pools.
    lhs_scale_mode : {"auto", "log", "linear"}
        Scale mapping for LHS candidates (passed through to LHS builder).
    xgb_use_gpu : bool or None
        GPU flag forwarded to XGBoost candidate builders.
    cuml_use_gpu : bool or None
        GPU flag forwarded to cuML SVM candidate builders.
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
            cuml_use_gpu=cuml_use_gpu,
            scale_mode=lhs_scale_mode,
            verbose=verbose,
            model_n_jobs=model_n_jobs,
        )
    return choose_param_grid(
        model_name=model_name,
        task=task,
        n_samples=n_samples,
        n_features=n_features,
        grid_size=budget,
        xgb_use_gpu=xgb_use_gpu,
        cuml_use_gpu=cuml_use_gpu,
        verbose=verbose,
        model_n_jobs=model_n_jobs,
    )


def _use_lhs_strategy(search_strategy: str, model_name: str) -> bool:
    """Return True when LHS sampling should be used instead of a fixed grid."""
    if search_strategy == "random":
        return True
    if search_strategy == "grid":
        return False
    # Hybrid: use LHS for high-dimensional tree/kernel/neural-network spaces.
    return model_name in {"rf", "xgb", "svm", "mlp"}


# ---------------------------------------------------------------------------
# Fixed-grid implementation
# ---------------------------------------------------------------------------


def choose_param_grid(
    model_name: str,
    task: str,
    n_samples: int,
    n_features: int,
    grid_size: int,
    xgb_use_gpu: bool | None = None,
    cuml_use_gpu: bool | None = None,
    verbose: bool = False,
    model_n_jobs: int = 1,
) -> list[dict[str, Any]]:
    """Build heuristic fixed-grid parameter candidates for GridSearchCV.

    Each returned dict maps pipeline parameter names to lists of candidate
    values. The "model" key holds the estimator instance(s) to try.
    """
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.linear_model import LinearRegression, LogisticRegression

    from .xgb_utils import xgb_gpu_available

    small = n_samples < 5000
    vlog(
        verbose,
        f"Choosing param grid for model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}, grid_size={grid_size}",
    )

    if model_name == "linear" and task == "classification":
        return [
            {
                "model": [
                    LogisticRegression(
                        C=np.inf,
                        solver="lbfgs",
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=RNG_SEED,
                    )
                ],
            }
        ]

    if model_name == "linear" and task == "regression":
        return [{"model": [LinearRegression(n_jobs=1)]}]

    if model_name == "rf":
        depth_base = max(3, int(np.log2(max(n_samples, 2))))
        axis_count = _axis_count_from_budget(grid_size, n_dims=4)
        n_estimators_pool = [80, 120, 180, 250] if small else [120, 200, 320, 480]
        depth_pool = [
            max(2, depth_base - 4),
            max(2, depth_base - 3),
            max(2, depth_base - 2),
            max(2, depth_base - 1),
            depth_base,
        ]
        min_leaf_pool = [4, 8, 16, 32]
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
        depth_low, depth_high = (2, 6) if small else (2, 8)
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

    if model_name == "svm":
        from .cuml_utils import cuml_gpu_available, get_cuml_svm

        use_gpu = bool(cuml_use_gpu) if cuml_use_gpu is not None else cuml_gpu_available()
        if use_gpu:
            base_model = get_cuml_svm(task)
        else:
            from sklearn.svm import LinearSVC, LinearSVR

            base_model = (
                LinearSVC(dual="auto", class_weight="balanced", max_iter=5000)
                if task == "classification"
                else LinearSVR(dual="auto", epsilon=0.0, max_iter=5000)
            )
        _C_SVM = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
        axis_count = _axis_count_from_budget(grid_size, n_dims=1)
        c_vals = _pick_evenly_spaced(_C_SVM, max(3, axis_count))
        return [{"model": [base_model], "model__C": c_vals}]

    if model_name == "mlp":
        from .deep import MLPClassifier, MLPRegressor

        _HIDDEN_PRESETS = [
            (64,),
            (128,),
            (128, 64),
            (256, 128),
            (256, 128, 64),
        ]
        _DROPOUT_POOL = [0.1, 0.2, 0.3, 0.4]
        _LR_POOL = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
        _WD_POOL = [1e-5, 1e-4, 1e-3]
        _BS_POOL = [64, 128, 256, 512]

        axis_count = _axis_count_from_budget(grid_size, n_dims=5)
        hidden_choices = _pick_evenly_spaced(_HIDDEN_PRESETS, axis_count)
        dropout_vals = _pick_evenly_spaced(_DROPOUT_POOL, axis_count)
        lr_vals = _pick_evenly_spaced(_LR_POOL, axis_count)
        wd_vals = _pick_evenly_spaced(_WD_POOL, axis_count)
        bs_vals = _pick_evenly_spaced(_BS_POOL, axis_count)

        base_model = (
            MLPClassifier(random_state=RNG_SEED)
            if task == "classification"
            else MLPRegressor(random_state=RNG_SEED)
        )
        return [
            {
                "model": [base_model],
                "model__hidden_sizes": hidden_choices,
                "model__dropout_rate": dropout_vals,
                "model__learning_rate": lr_vals,
                "model__weight_decay": wd_vals,
                "model__batch_size": bs_vals,
            }
        ]

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
    cuml_use_gpu: bool | None = None,
    scale_mode: str = "auto",
    verbose: bool = False,
    model_n_jobs: int = 1,
) -> list[dict[str, list[Any]]]:
    """Build a space-filling LHS candidate set for GridSearchCV.

    Returns a list of singleton dicts compatible with
    ``GridSearchCV(param_grid=<result>)``.
    """
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.linear_model import LinearRegression, LogisticRegression

    from .xgb_utils import xgb_gpu_available

    n_points = max(1, int(budget))
    small = n_samples < 5000
    vlog(
        verbose,
        f"Choosing LHS candidate space model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}, "
        f"budget={n_points}, scale_mode={scale_mode}",
    )

    if model_name == "linear" and task == "classification":
        return [
            {
                "model": [
                    LogisticRegression(
                        C=np.inf,
                        solver="lbfgs",
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=RNG_SEED,
                    )
                ],
            }
        ]

    if model_name == "linear" and task == "regression":
        return [{"model": [LinearRegression(n_jobs=1)]}]

    if model_name == "rf":
        depth_base = max(3, int(np.log2(max(n_samples, 2))))
        depth_low, depth_high = max(2, depth_base - 4), depth_base
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
                                    4.0,
                                    32.0,
                                    scale_mode,
                                    default_scale="log",
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
        depth_low, depth_high = (2, 6) if small else (2, 8)
        n_est_low, n_est_high = (60, 200) if small else (100, 500)
        unit = _lhs_unit(n_points=n_points, n_dims=6, seed=RNG_SEED)
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
                    "model__min_child_weight": [
                        int(
                            round(
                                _map_with_scale(
                                    float(row[5]), 1.0, 20.0, scale_mode, default_scale="log"
                                )
                            )
                        )
                    ],
                }
            )
        return out

    if model_name == "svm":
        from .cuml_utils import cuml_gpu_available, get_cuml_svm

        use_gpu = bool(cuml_use_gpu) if cuml_use_gpu is not None else cuml_gpu_available()
        if use_gpu:
            base_model = get_cuml_svm(task)
        else:
            from sklearn.svm import LinearSVC, LinearSVR

            base_model = (
                LinearSVC(dual="auto", class_weight="balanced", max_iter=5000)
                if task == "classification"
                else LinearSVR(dual="auto", epsilon=0.0, max_iter=5000)
            )
        # 1 LHS dim: C on log scale — same range for both tasks
        unit = _lhs_unit(n_points=n_points, n_dims=1, seed=RNG_SEED)
        return [
            {
                "model": [base_model],
                "model__C": [
                    _map_with_scale(float(row[0]), 0.001, 100.0, scale_mode, default_scale="log")
                ],
            }
            for row in unit
        ]

    if model_name == "mlp":
        from .deep import MLPClassifier, MLPRegressor

        # 5D LHS space:
        #   dim 0 — hidden_sizes  (index into preset architectures, linear)
        #   dim 1 — dropout_rate  [0.05, 0.50], linear
        #   dim 2 — learning_rate [1e-4, 1e-2], log
        #   dim 3 — weight_decay  [1e-6, 1e-3], log
        #   dim 4 — batch_size    (index into preset sizes, linear)
        _HIDDEN_PRESETS = [
            (64,),
            (128,),
            (128, 64),
            (256, 128),
            (256, 128, 64),
        ]
        _BATCH_SIZES = [64, 128, 256, 512]

        base_model = (
            MLPClassifier(random_state=RNG_SEED)
            if task == "classification"
            else MLPRegressor(random_state=RNG_SEED)
        )
        unit = _lhs_unit(n_points=n_points, n_dims=5, seed=RNG_SEED)
        out: list[dict[str, list[Any]]] = []
        for row in unit:
            hidden_idx = min(int(float(row[0]) * len(_HIDDEN_PRESETS)), len(_HIDDEN_PRESETS) - 1)
            bs_idx = min(int(float(row[4]) * len(_BATCH_SIZES)), len(_BATCH_SIZES) - 1)
            out.append(
                {
                    "model": [base_model],
                    "model__hidden_sizes": [_HIDDEN_PRESETS[hidden_idx]],
                    "model__dropout_rate": [
                        _map_with_scale(
                            float(row[1]), 0.05, 0.50, scale_mode, default_scale="linear"
                        )
                    ],
                    "model__learning_rate": [
                        _map_with_scale(
                            float(row[2]), 1e-4, 1e-2, scale_mode, default_scale="log"
                        )
                    ],
                    "model__weight_decay": [
                        _map_with_scale(
                            float(row[3]), 1e-6, 1e-3, scale_mode, default_scale="log"
                        )
                    ],
                    "model__batch_size": [_BATCH_SIZES[bs_idx]],
                }
            )
        return out

    raise ValueError(f"Unsupported model name: {model_name}")
