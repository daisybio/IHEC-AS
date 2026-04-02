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

import os

import numpy as np

from ..config import RNG_SEED
from ..utils import vlog

__all__ = [
    "build_param_candidates",
    "build_optuna_distributions",
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
#
# Grid sizing tradeoff: LogisticRegressionCV parallelises over
# n_l1_ratios × inner_folds tasks.  Reducing l1_ratio count cuts parallel
# batches (structural); reducing alpha/C count halves per-task work (warm-start
# path means all Cs are solved in one pass per l1_ratio per fold).
#
# Reduced grids (default): 50 alpha/C values, 5 l1_ratios → ~3 batches at 16c
# Original grids (commented): 100 alpha/C values, 7 l1_ratios → ~4 batches at 16c
_ALPHA_GRID: np.ndarray = np.logspace(-6, 2, num=50)
_C_GRID: np.ndarray = 1.0 / _ALPHA_GRID
_L1_RATIO_GRID: list[float] = [0.1, 0.5, 0.9, 0.95, 1.0]
# UNUSED: Denser legacy grid retained as a reference for historical tuning behavior.
# TODO: remove this block if benchmark history is captured elsewhere.
# Original (denser) grids — ~2× slower per fold, marginally better coverage:
# _ALPHA_GRID: np.ndarray = np.logspace(-6, 2, num=100)
# _C_GRID: np.ndarray = 1.0 / _ALPHA_GRID
# _L1_RATIO_GRID: list[float] = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]


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
    scale_pos_weight: float | None = None,
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
            scale_pos_weight=scale_pos_weight,
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
        scale_pos_weight=scale_pos_weight,
    )


def _use_lhs_strategy(search_strategy: str, model_name: str) -> bool:
    """Return True when LHS sampling should be used instead of a fixed grid."""
    if search_strategy == "random":
        return True
    if search_strategy == "grid":
        return False
    # Hybrid: use LHS for high-dimensional tree/kernel/neural-network spaces.
    return model_name in {"rf", "xgb", "svm", "nysvm", "mlp"}


def _resolve_gpu_flag(requested: bool | None, detect_fn: Any) -> bool:
    """Resolve explicit/auto GPU flag.

    If ``requested`` is not None, respect it; otherwise call ``detect_fn``.
    """
    return bool(requested) if requested is not None else bool(detect_fn())


def _linear_candidates(
    task: str, cuml_use_gpu: bool | None
) -> list[dict[str, list[Any]]]:
    """Return single-candidate model dicts for linear models."""
    from sklearn.linear_model import LinearRegression, LogisticRegression

    from .cuml_utils import cuml_gpu_available, get_cuml_linear

    use_gpu = _resolve_gpu_flag(cuml_use_gpu, cuml_gpu_available)
    if use_gpu:
        return [{"model": [get_cuml_linear(task)]}]
    if task == "classification":
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
    return [{"model": [LinearRegression(n_jobs=1)]}]


def _svm_base_model(task: str, cuml_use_gpu: bool | None) -> Any:
    """Return SVM estimator (cuML if available/requested, else sklearn)."""
    from .cuml_utils import cuml_gpu_available, get_cuml_svm

    use_gpu = _resolve_gpu_flag(cuml_use_gpu, cuml_gpu_available)
    if use_gpu:
        return get_cuml_svm(task)

    from sklearn.svm import SVC, SVR

    return (
        SVC(
            kernel="rbf",
            class_weight="balanced",
            probability=True,
            max_iter=5000,
            tol=1e-3,
            cache_size=2000,
        )
        if task == "classification"
        else SVR(kernel="rbf", max_iter=5000, tol=1e-3, cache_size=2000)
    )


def _xgb_common_kwargs(use_gpu: bool, model_n_jobs: int) -> dict[str, Any]:
    """Common kwargs shared by XGBoost estimators."""
    return {
        "random_state": RNG_SEED,
        "n_jobs": 1 if use_gpu else model_n_jobs,
        "tree_method": "hist",
        "device": "cuda" if use_gpu else "cpu",
    }


def _xgbrf_common_kwargs(use_gpu: bool, model_n_jobs: int) -> dict[str, Any]:
    """Common kwargs for XGBRF estimators."""
    common = _xgb_common_kwargs(use_gpu, model_n_jobs)
    common["subsample"] = 0.8
    return common


def _build_xgbrf_base_model(
    task: str,
    common: dict[str, Any],
    scale_pos_weight: float | None,
) -> Any:
    """Construct XGBRF classifier/regressor base estimator."""
    try:
        from xgboost import XGBRFClassifier, XGBRFRegressor
    except Exception as exc:
        raise RuntimeError("xgboost requested but not installed") from exc

    if task == "classification":
        return XGBRFClassifier(
            eval_metric="logloss",
            scale_pos_weight=scale_pos_weight,
            **common,
        )
    return XGBRFRegressor(eval_metric="rmse", **common)


def _build_xgb_base_model(task: str, common: dict[str, Any]) -> Any:
    """Construct XGBoost classifier/regressor base estimator."""
    try:
        from xgboost import XGBClassifier, XGBRegressor
    except Exception as exc:
        raise RuntimeError("xgboost requested but not installed") from exc

    return (
        XGBClassifier(eval_metric="logloss", **common)
        if task == "classification"
        else XGBRegressor(eval_metric="rmse", **common)
    )


def _rf_ranges(n_samples: int) -> dict[str, float]:
    """Shared RF/XGBRF heuristic ranges from dataset size."""
    log10_n = max(1.0, float(np.log10(max(n_samples, 2))))
    depth_cap = max(8, min(12, int(36 - 5.5 * log10_n)))
    depth_base = min(max(3, int(np.log2(max(n_samples, 2)))), depth_cap)
    n_est_high = max(150, int(300 - 30 * max(0.0, log10_n - 3.5)))
    min_child_low = max(1, int(log10_n - 2) + 1)
    min_child_high = max(10, int(log10_n**2 - 5))
    return {
        "log10_n": log10_n,
        "depth_base": float(depth_base),
        "n_est_low": 100.0,
        "n_est_high": float(n_est_high),
        "min_child_low": float(min_child_low),
        "min_child_high": float(min_child_high),
    }


def _xgb_ranges(n_samples: int) -> dict[str, float]:
    """Shared XGBoost heuristic ranges from dataset size."""
    log2_n = max(1.0, float(np.log2(max(n_samples, 2))))
    log10_n = max(1.0, float(np.log10(max(n_samples, 2))))
    depth_cap = min(7, max(3, int(log2_n // 3)))
    depth_low = max(2, depth_cap - 2)
    depth_high = depth_cap
    lr_low = round(min(0.05, max(0.01, 0.005 * log10_n)), 3)
    lr_high = round(min(0.15, max(0.06, 0.65 / log10_n)), 3)
    n_est_high = max(150, int(300 - 30 * max(0.0, log10_n - 3.5)))
    min_child_low = max(1, int(log10_n - 2) + 1)
    min_child_high = max(10, int(log10_n**2 - 5))
    return {
        "depth_low": float(depth_low),
        "depth_high": float(depth_high),
        "lr_low": float(lr_low),
        "lr_high": float(lr_high),
        "n_est_low": 100.0,
        "n_est_high": float(n_est_high),
        "min_child_low": float(min_child_low),
        "min_child_high": float(min_child_high),
        "subsample_low": 0.7,
        "subsample_high": 0.95,
        "colsample_low": 0.5,
        "colsample_high": 0.9,
    }


def _svm_ranges(n_samples: int) -> tuple[float, float, float, float]:
    """Return (c_low, c_high, gamma_low, gamma_high) for SVM."""
    c_low = 0.03
    c_high = round(max(10.0, min(100.0, 3000.0 * float(n_samples) ** (-0.4))), 1)
    gamma_low = 0.0003
    gamma_high = round(max(0.001, min(0.1, 3.0 * float(n_samples) ** (-0.4))), 4)
    return c_low, c_high, gamma_low, gamma_high


def _available_memory_bytes() -> int | None:
    """Best-effort available system memory in bytes."""
    try:
        import psutil

        return int(psutil.virtual_memory().available)
    except Exception:
        pass

    # Linux fallback for HPC nodes where psutil may be unavailable.
    try:
        with open("/proc/meminfo", "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return int(parts[1]) * 1024
    except Exception:
        pass

    try:
        pages = int(os.sysconf("SC_AVPHYS_PAGES"))
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        return pages * page_size
    except Exception:
        return None


def _nysvm_ranges(n_samples: int) -> dict[str, Any]:
    """Compute Nyström SVM parameter bounds from sample count.

    n_components controls approximation quality at O(n × k) memory cost per
    trial.  We combine sample-count scaling with a memory-aware cap derived
    from currently available RAM to avoid OOM kills on large folds.

    gamma uses the same adaptive formula as ``_svm_ranges`` — the RBF
    bandwidth interpretation is identical for Nyström and full kernel SVM.
    A fixed gamma_high of 1.0 would be 40–100× above the useful range for
    these dataset sizes and produce a near-identity kernel.
    """
    # Nyström requires n_components <= n_samples; keep ranges fold-safe.
    n_cap = max(1, int(n_samples))
    n_comp_low = min(max(32, min(200, n_samples // 3000)), n_cap)
    n_comp_high = min(max(96, min(500, n_samples // 2000)), n_cap)

    # Optional manual cap for cluster tuning, e.g. SPLICING_ML_NYSVM_MAX_COMPONENTS=192.
    env_cap_raw = os.getenv("SPLICING_ML_NYSVM_MAX_COMPONENTS")
    if env_cap_raw:
        try:
            n_comp_high = min(n_comp_high, max(32, int(env_cap_raw)))
        except ValueError:
            pass

    avail_bytes = _available_memory_bytes()
    if avail_bytes is not None and n_samples > 0:
        # Approximation for Nyström feature map memory in float64 plus workspace.
        # Reserve only a fraction of available RAM to account for CV/model overhead.
        bytes_per_value = 8.0
        workspace_factor = 2.0
        usable_fraction = 0.12
        mem_cap = int(
            (float(avail_bytes) * usable_fraction)
            / (float(n_samples) * bytes_per_value * workspace_factor)
        )
        n_comp_high = min(n_comp_high, max(32, mem_cap))

    if n_comp_low > n_comp_high:
        n_comp_low = n_comp_high
    gamma_low = 0.0003
    gamma_high = round(max(0.003, min(0.3, 3.0 * float(n_samples) ** (-0.4))), 4)
    return {
        "n_comp_low": n_comp_low,
        "n_comp_high": n_comp_high,
        "gamma_low": gamma_low,
        "gamma_high": gamma_high,
    }


def _rf_candidate_from_row(
    base_model: Any,
    row: np.ndarray,
    *,
    n_est_low: int,
    n_est_high: int,
    depth_low: int,
    depth_high: int,
    min_child_low: float,
    min_child_high: float,
    scale_mode: str,
) -> dict[str, list[Any]]:
    """Build one RF/XGBRF candidate dict from a 4D sampled row."""
    return {
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
        "model__colsample_bynode": [
            _map_with_scale(
                float(row[2]),
                0.3,
                0.8,
                scale_mode,
                default_scale="linear",
            )
        ],
        "model__min_child_weight": [
            int(
                round(
                    _map_with_scale(
                        float(row[3]),
                        min_child_low,
                        min_child_high,
                        scale_mode,
                        default_scale="log",
                    )
                )
            )
        ],
    }


def _xgb_candidate_from_row(
    base_model: Any,
    row: np.ndarray,
    *,
    n_est_low: int,
    n_est_high: int,
    depth_low: int,
    depth_high: int,
    lr_low: float,
    lr_high: float,
    subsample_low: float,
    subsample_high: float,
    colsample_low: float,
    colsample_high: float,
    scale_mode: str,
    min_child_low: float | None = None,
    min_child_high: float | None = None,
) -> dict[str, list[Any]]:
    """Build one XGB candidate dict from a sampled row.

    Uses the first 5 dims for base XGB params. If both min-child bounds are
    provided, dim 5 is used for ``model__min_child_weight``.
    """
    out: dict[str, list[Any]] = {
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
                float(row[2]),
                lr_low,
                lr_high,
                scale_mode,
                default_scale="log",
            )
        ],
        "model__subsample": [
            _map_with_scale(
                float(row[3]),
                subsample_low,
                subsample_high,
                scale_mode,
                default_scale="linear",
            )
        ],
        "model__colsample_bytree": [
            _map_with_scale(
                float(row[4]),
                colsample_low,
                colsample_high,
                scale_mode,
                default_scale="linear",
            )
        ],
    }
    if min_child_low is not None and min_child_high is not None:
        out["model__min_child_weight"] = [
            int(
                round(
                    _map_with_scale(
                        float(row[5]),
                        min_child_low,
                        min_child_high,
                        scale_mode,
                        default_scale="log",
                    )
                )
            )
        ]
    return out


def _mlp_candidate_from_row(
    base_model: Any,
    row: np.ndarray,
    scale_mode: str,
    batch_size: int,
) -> dict[str, list[Any]]:
    """Build one MLP candidate dict from a 4D sampled row.

    Dimensions: [hidden_sizes (linear), dropout_rate (linear), learning_rate (log), weight_decay (log)].
    Batch size is fixed, not tuned.
    """
    hidden_idx = min(
        int(float(row[0]) * len(_MLP_HIDDEN_PRESETS)),
        len(_MLP_HIDDEN_PRESETS) - 1,
    )
    return {
        "model": [base_model],
        "model__hidden_sizes": [_MLP_HIDDEN_PRESETS[hidden_idx]],
        "model__dropout_rate": [
            _map_with_scale(
                float(row[1]),
                0.20,
                0.55,
                scale_mode,
                default_scale="linear",
            )
        ],
        "model__learning_rate": [
            _map_with_scale(
                float(row[2]),
                1e-4,
                3e-3,
                scale_mode,
                default_scale="log",
            )
        ],
        "model__weight_decay": [
            _map_with_scale(
                float(row[3]),
                1e-5,
                5e-3,
                scale_mode,
                default_scale="log",
            )
        ],
        "model__batch_size": [batch_size],
    }


_MLP_HIDDEN_PRESETS = [
    (64,),
    (128,),
    (128, 64),
    (256, 128),
]
_MLP_DROPOUT_POOL = [0.2, 0.3, 0.4, 0.5]
_MLP_LR_POOL = [1e-4, 3e-4, 1e-3, 3e-3]
_MLP_WD_POOL = [1e-5, 1e-4, 1e-3, 3e-3]
_MLP_BATCH_SIZES = [128, 256, 512]


def _mlp_batch_sizes(n_samples: int) -> list[int]:
    """Return batch size pool scaled to n_samples.

    Larger datasets benefit from larger batches: fewer gradient steps per
    epoch dramatically reduces wall time with negligible accuracy impact at
    this scale.  Rule of thumb: target ~100–300 batches/epoch.
    """
    if n_samples < 20_000:
        return [128, 256, 512]
    elif n_samples < 100_000:
        return [256, 512, 1024]
    elif n_samples < 400_000:
        return [512, 1024, 2048]
    else:
        return [1024, 2048, 4096]


def _mlp_base_model(task: str) -> Any:
    """Return MLP estimator for the task."""
    from .deep import MLPClassifier, MLPRegressor

    return (
        MLPClassifier(random_state=RNG_SEED)
        if task == "classification"
        else MLPRegressor(random_state=RNG_SEED)
    )


def _nysvm_base_model(task: str, n_components: int = 300) -> Any:
    """Return Nyström-approximated linear SVM as a nested sklearn Pipeline.

    The pipeline has two steps:
      - ``nystroem``: maps input features to an approximate RBF kernel space
      - ``svc`` / ``svr``: LinearSVC/SVR fitted in that feature space (O(n))

    Parameter keys for grid/Optuna search:
      ``model__nystroem__n_components``, ``model__nystroem__gamma``,
      ``model__svc__C`` (classification) or ``model__svr__C`` (regression).
    """
    from sklearn.kernel_approximation import Nystroem
    from sklearn.pipeline import Pipeline as _SkPipeline
    from sklearn.svm import LinearSVC, LinearSVR

    nystroem = Nystroem(kernel="rbf", n_components=n_components, random_state=RNG_SEED)
    if task == "classification":
        return _SkPipeline(
            [
                ("nystroem", nystroem),
                ("svc", LinearSVC(class_weight="balanced", max_iter=20000)),
            ]
        )
    return _SkPipeline([("nystroem", nystroem), ("svr", LinearSVR(max_iter=20000))])


def _nysvm_candidate_from_row(
    base_model: Any,
    row: np.ndarray,
    *,
    n_comp_low: int,
    n_comp_high: int,
    gamma_low: float,
    gamma_high: float,
    c_key: str,
    scale_mode: str,
) -> dict[str, list[Any]]:
    """Build one Nyström SVM candidate dict from a 3D LHS row.

    Dimensions: [n_components (linear), gamma (log), C (log)].
    ``c_key`` is ``"model__svc__C"`` or ``"model__svr__C"`` depending on task.
    """
    return {
        "model": [base_model],
        "model__nystroem__n_components": [
            int(
                round(
                    _map_with_scale(
                        float(row[0]),
                        n_comp_low,
                        n_comp_high,
                        scale_mode,
                        default_scale="linear",
                    )
                )
            )
        ],
        "model__nystroem__gamma": [
            _map_with_scale(
                float(row[1]), gamma_low, gamma_high, scale_mode, default_scale="log"
            )
        ],
        c_key: [
            _map_with_scale(float(row[2]), 0.01, 100.0, scale_mode, default_scale="log")
        ],
    }


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
    scale_pos_weight: float | None = None,
) -> list[dict[str, Any]]:
    """Build heuristic fixed-grid parameter candidates for GridSearchCV.

    Each returned dict maps pipeline parameter names to lists of candidate
    values. The "model" key holds the estimator instance(s) to try.
    """
    from .xgb_utils import xgb_gpu_available

    vlog(
        verbose,
        f"Choosing param grid for model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}, grid_size={grid_size}",
    )

    if model_name == "linear":
        return _linear_candidates(task, cuml_use_gpu)

    if model_name == "rf":
        use_gpu = _resolve_gpu_flag(xgb_use_gpu, xgb_gpu_available)
        common = _xgbrf_common_kwargs(use_gpu, model_n_jobs)
        ranges = _rf_ranges(n_samples)
        depth_base = int(ranges["depth_base"])
        n_est_high = int(ranges["n_est_high"])
        min_child_low = int(ranges["min_child_low"])
        min_child_high = int(ranges["min_child_high"])
        axis_count = _axis_count_from_budget(grid_size, n_dims=4)
        n_estimators_pool = [80, 100, int((100 + n_est_high) / 2), n_est_high]
        depth_pool = [max(2, depth_base - i) for i in range(4, -1, -1)]
        colsample_pool = [0.3, 0.5, 0.7]
        min_child_pool = [
            min_child_low,
            max(min_child_low + 1, int(min_child_low * 3)),
            max(min_child_low + 2, int(min_child_high * 0.5)),
            min_child_high,
        ]
        n_estimators = _pick_evenly_spaced(n_estimators_pool, axis_count)
        max_depth = _pick_evenly_spaced(depth_pool, axis_count)
        colsample = _pick_evenly_spaced(colsample_pool, axis_count)
        min_child = _pick_evenly_spaced(min_child_pool, axis_count)
        base_model = _build_xgbrf_base_model(task, common, scale_pos_weight)
        return [
            {
                "model": [base_model],
                "model__n_estimators": n_estimators,
                "model__max_depth": max_depth,
                "model__colsample_bynode": colsample,
                "model__min_child_weight": min_child,
            }
        ]

    if model_name == "xgb":
        use_gpu = _resolve_gpu_flag(xgb_use_gpu, xgb_gpu_available)
        common = _xgb_common_kwargs(use_gpu, model_n_jobs)
        n_points = max(1, int(grid_size))
        ranges = _xgb_ranges(n_samples)
        depth_low = int(ranges["depth_low"])
        depth_high = int(ranges["depth_high"])
        lr_low = ranges["lr_low"]
        lr_high = ranges["lr_high"]
        n_est_low = int(ranges["n_est_low"])
        n_est_high = int(ranges["n_est_high"])
        subsample_low = ranges["subsample_low"]
        subsample_high = ranges["subsample_high"]
        colsample_low = ranges["colsample_low"]
        colsample_high = ranges["colsample_high"]
        unit = _lhs_unit(n_points=n_points, n_dims=5, seed=RNG_SEED)
        base_model = _build_xgb_base_model(task, common)
        return [
            _xgb_candidate_from_row(
                base_model,
                row,
                n_est_low=n_est_low,
                n_est_high=n_est_high,
                depth_low=depth_low,
                depth_high=depth_high,
                lr_low=lr_low,
                lr_high=lr_high,
                subsample_low=subsample_low,
                subsample_high=subsample_high,
                colsample_low=colsample_low,
                colsample_high=colsample_high,
                scale_mode="linear",
            )
            for row in unit
        ]

    if model_name == "svm":
        base_model = _svm_base_model(task, cuml_use_gpu)
        _, c_high, _, gamma_high = _svm_ranges(n_samples)
        _C_SVM = [
            v for v in [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0] if v <= c_high
        ]
        _GAMMA_SVM = [
            v for v in [0.0003, 0.001, 0.003, 0.01, 0.03, 0.1] if v <= gamma_high
        ]
        axis_count = _axis_count_from_budget(grid_size, n_dims=2)
        c_vals = _pick_evenly_spaced(_C_SVM, max(3, axis_count))
        gamma_vals = _pick_evenly_spaced(_GAMMA_SVM, max(2, axis_count))
        return [{"model": [base_model], "model__C": c_vals, "model__gamma": gamma_vals}]

    if model_name == "nysvm":
        ranges = _nysvm_ranges(n_samples)
        n_comp_low, n_comp_high = ranges["n_comp_low"], ranges["n_comp_high"]
        gamma_high = ranges["gamma_high"]
        base_model = _nysvm_base_model(task, (n_comp_low + n_comp_high) // 2)
        c_key = "model__svc__C" if task == "classification" else "model__svr__C"
        step = max(1, (n_comp_high - n_comp_low) // 3)
        n_comp_pool = sorted(
            {n_comp_low, n_comp_low + step, n_comp_high - step, n_comp_high}
        )
        gamma_pool = [
            v for v in [0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3] if v <= gamma_high
        ]
        axis_count = _axis_count_from_budget(grid_size, n_dims=3)
        return [
            {
                "model": [base_model],
                "model__nystroem__n_components": _pick_evenly_spaced(
                    n_comp_pool, axis_count
                ),
                "model__nystroem__gamma": _pick_evenly_spaced(gamma_pool, axis_count),
                c_key: _pick_evenly_spaced(
                    [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0], axis_count
                ),
            }
        ]

    if model_name == "mlp":
        axis_count = _axis_count_from_budget(grid_size, n_dims=4)
        hidden_choices = _pick_evenly_spaced(_MLP_HIDDEN_PRESETS, axis_count)
        dropout_vals = _pick_evenly_spaced(_MLP_DROPOUT_POOL, axis_count)
        lr_vals = _pick_evenly_spaced(_MLP_LR_POOL, axis_count)
        wd_vals = _pick_evenly_spaced(_MLP_WD_POOL, axis_count)
        batch_sizes = _mlp_batch_sizes(n_samples)
        batch_size = batch_sizes[len(batch_sizes) // 2]

        base_model = _mlp_base_model(task)
        return [
            {
                "model": [base_model],
                "model__hidden_sizes": hidden_choices,
                "model__dropout_rate": dropout_vals,
                "model__learning_rate": lr_vals,
                "model__weight_decay": wd_vals,
                "model__batch_size": [batch_size],
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
    scale_pos_weight: float | None = None,
) -> list[dict[str, list[Any]]]:
    """Build a space-filling LHS candidate set for GridSearchCV.

    Returns a list of singleton dicts compatible with
    ``GridSearchCV(param_grid=<result>)``.
    """
    from .xgb_utils import xgb_gpu_available

    n_points = max(1, int(budget))
    vlog(
        verbose,
        f"Choosing LHS candidate space model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}, "
        f"budget={n_points}, scale_mode={scale_mode}",
    )

    if model_name == "linear":
        return _linear_candidates(task, cuml_use_gpu)

    if model_name == "rf":
        use_gpu = _resolve_gpu_flag(xgb_use_gpu, xgb_gpu_available)
        common = _xgbrf_common_kwargs(use_gpu, model_n_jobs)
        ranges = _rf_ranges(n_samples)
        depth_base = int(ranges["depth_base"])
        depth_low = max(2, depth_base - 4)
        depth_high = depth_base
        n_est_low = int(ranges["n_est_low"])
        n_est_high = int(ranges["n_est_high"])
        min_child_low = ranges["min_child_low"]
        min_child_high = ranges["min_child_high"]
        base_model = _build_xgbrf_base_model(task, common, scale_pos_weight)
        unit = _lhs_unit(n_points=n_points, n_dims=4, seed=RNG_SEED)
        return [
            _rf_candidate_from_row(
                base_model,
                row,
                n_est_low=n_est_low,
                n_est_high=n_est_high,
                depth_low=depth_low,
                depth_high=depth_high,
                min_child_low=min_child_low,
                min_child_high=min_child_high,
                scale_mode=scale_mode,
            )
            for row in unit
        ]

    if model_name == "xgb":
        use_gpu = _resolve_gpu_flag(xgb_use_gpu, xgb_gpu_available)
        common = _xgb_common_kwargs(use_gpu, model_n_jobs)
        base_model = _build_xgb_base_model(task, common)
        ranges = _xgb_ranges(n_samples)
        depth_low = int(ranges["depth_low"])
        depth_high = int(ranges["depth_high"])
        lr_low = ranges["lr_low"]
        lr_high = ranges["lr_high"]
        n_est_low = int(ranges["n_est_low"])
        n_est_high = int(ranges["n_est_high"])
        min_child_low = ranges["min_child_low"]
        min_child_high = ranges["min_child_high"]
        subsample_low = ranges["subsample_low"]
        subsample_high = ranges["subsample_high"]
        colsample_low = ranges["colsample_low"]
        colsample_high = ranges["colsample_high"]
        unit = _lhs_unit(n_points=n_points, n_dims=6, seed=RNG_SEED)
        return [
            _xgb_candidate_from_row(
                base_model,
                row,
                n_est_low=n_est_low,
                n_est_high=n_est_high,
                depth_low=depth_low,
                depth_high=depth_high,
                lr_low=lr_low,
                lr_high=lr_high,
                subsample_low=subsample_low,
                subsample_high=subsample_high,
                colsample_low=colsample_low,
                colsample_high=colsample_high,
                scale_mode=scale_mode,
                min_child_low=min_child_low,
                min_child_high=min_child_high,
            )
            for row in unit
        ]

    if model_name == "svm":
        base_model = _svm_base_model(task, cuml_use_gpu)
        # 2 LHS dims: C and gamma on log scale
        c_low, c_high, gamma_low, gamma_high = _svm_ranges(n_samples)
        unit = _lhs_unit(n_points=n_points, n_dims=2, seed=RNG_SEED)
        return [
            {
                "model": [base_model],
                "model__C": [
                    _map_with_scale(
                        float(row[0]), c_low, c_high, scale_mode, default_scale="log"
                    )
                ],
                "model__gamma": [
                    _map_with_scale(
                        float(row[1]),
                        gamma_low,
                        gamma_high,
                        scale_mode,
                        default_scale="log",
                    )
                ],
            }
            for row in unit
        ]

    if model_name == "nysvm":
        # 3D LHS space: n_components (linear), gamma (log), C (log)
        ranges = _nysvm_ranges(n_samples)
        n_comp_low = ranges["n_comp_low"]
        n_comp_high = ranges["n_comp_high"]
        gamma_low = ranges["gamma_low"]
        gamma_high = ranges["gamma_high"]
        base_model = _nysvm_base_model(task, (n_comp_low + n_comp_high) // 2)
        c_key = "model__svc__C" if task == "classification" else "model__svr__C"
        unit = _lhs_unit(n_points=n_points, n_dims=3, seed=RNG_SEED)
        return [
            _nysvm_candidate_from_row(
                base_model,
                row,
                n_comp_low=n_comp_low,
                n_comp_high=n_comp_high,
                gamma_low=gamma_low,
                gamma_high=gamma_high,
                c_key=c_key,
                scale_mode=scale_mode,
            )
            for row in unit
        ]

    if model_name == "mlp":
        # 4D LHS space (batch_size fixed to mid-range value, not tuned):
        #   dim 0 — hidden_sizes  (index into preset architectures, linear)
        #   dim 1 — dropout_rate  [0.20, 0.55], linear
        #   dim 2 — learning_rate [1e-4, 3e-3], log
        #   dim 3 — weight_decay  [1e-5, 5e-3], log
        base_model = _mlp_base_model(task)
        batch_sizes = _mlp_batch_sizes(n_samples)
        batch_size = batch_sizes[len(batch_sizes) // 2]  # mid-range value
        unit = _lhs_unit(n_points=n_points, n_dims=4, seed=RNG_SEED)
        return [
            _mlp_candidate_from_row(base_model, row, scale_mode, batch_size=batch_size)
            for row in unit
        ]

    raise ValueError(f"Unsupported model name: {model_name}")


# ---------------------------------------------------------------------------
# Optuna distribution builder
# ---------------------------------------------------------------------------


def build_optuna_distributions(
    model_name: str,
    task: str,
    n_samples: int,
    n_features: int,
    xgb_use_gpu: bool | None = None,
    cuml_use_gpu: bool | None = None,
    verbose: bool = False,
    model_n_jobs: int = 1,
    scale_pos_weight: float | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Return (base_estimator, optuna_distributions) for OptunaSearchCV.

    Maps each model's parameter space to Optuna distribution objects using the
    same numeric ranges as the LHS candidate builders.  The base_estimator is a
    fully configured (unfitted) sklearn estimator to set on the pipeline "model"
    step.

    Parameters
    ----------
    model_name : str
        One of: "rf", "xgb", "svm", "mlp".  "linear" is not supported
        (no hyperparameters to tune).
    task : str
        "regression" or "classification".
    n_samples, n_features : int
        Dataset dimensions used to select size-conditional parameter ranges.
    xgb_use_gpu : bool or None
        GPU flag forwarded to XGBoost/RF estimator constructors.
    cuml_use_gpu : bool or None
        GPU flag forwarded to cuML SVM estimator constructor.
    verbose : bool
        Enable debug logging.
    model_n_jobs : int
        Thread budget for the base estimator itself (not the search).
    scale_pos_weight : float or None
        Class imbalance weight forwarded to XGB classification estimators.

    Returns
    -------
    base_estimator : estimator
        Model instance (unfitted) to set as the pipeline "model" step.
    distributions : dict[str, optuna.distributions.*]
        Pipeline-namespaced parameter names mapped to Optuna distributions.

    Raises
    ------
    ImportError
        If ``optuna`` is not installed.
    ValueError
        If ``model_name`` is not supported.
    """
    try:
        from optuna.distributions import (
            CategoricalDistribution,
            FloatDistribution,
            IntDistribution,
        )
    except ImportError as exc:
        raise ImportError(
            "optuna_backend=True requires optuna and optuna-integration[sklearn]. "
            "Install with: pip install 'optuna-integration[sklearn]>=3.4'"
        ) from exc

    from .xgb_utils import xgb_gpu_available

    vlog(
        verbose,
        f"Building Optuna distributions for model={model_name}, task={task}, "
        f"n_samples={n_samples}, n_features={n_features}",
    )

    if model_name == "rf":
        use_gpu = _resolve_gpu_flag(xgb_use_gpu, xgb_gpu_available)
        common = _xgbrf_common_kwargs(use_gpu, model_n_jobs)
        ranges = _rf_ranges(n_samples)
        depth_base = int(ranges["depth_base"])
        depth_low = max(2, depth_base - 4)
        depth_high = depth_base
        n_est_low = int(ranges["n_est_low"])
        n_est_high = int(ranges["n_est_high"])
        min_child_low = int(ranges["min_child_low"])
        min_child_high = int(ranges["min_child_high"])
        base_model = _build_xgbrf_base_model(task, common, scale_pos_weight)
        distributions: dict[str, Any] = {
            "model__n_estimators": IntDistribution(n_est_low, n_est_high),
            "model__max_depth": IntDistribution(depth_low, depth_high),
            "model__colsample_bynode": FloatDistribution(0.3, 0.8),
            "model__min_child_weight": IntDistribution(
                min_child_low, min_child_high, log=True
            ),
        }
        return base_model, distributions

    if model_name == "xgb":
        use_gpu = _resolve_gpu_flag(xgb_use_gpu, xgb_gpu_available)
        common = _xgb_common_kwargs(use_gpu, model_n_jobs)
        ranges = _xgb_ranges(n_samples)
        depth_low = int(ranges["depth_low"])
        depth_high = int(ranges["depth_high"])
        lr_low = ranges["lr_low"]
        lr_high = ranges["lr_high"]
        n_est_low = int(ranges["n_est_low"])
        n_est_high = int(ranges["n_est_high"])
        min_child_low = int(ranges["min_child_low"])
        min_child_high = int(ranges["min_child_high"])
        subsample_low = ranges["subsample_low"]
        subsample_high = ranges["subsample_high"]
        colsample_low = ranges["colsample_low"]
        colsample_high = ranges["colsample_high"]
        base_model = _build_xgb_base_model(task, common)
        distributions = {
            "model__n_estimators": IntDistribution(n_est_low, n_est_high),
            "model__max_depth": IntDistribution(depth_low, depth_high),
            "model__learning_rate": FloatDistribution(lr_low, lr_high, log=True),
            "model__subsample": FloatDistribution(subsample_low, subsample_high),
            "model__colsample_bytree": FloatDistribution(colsample_low, colsample_high),
            "model__min_child_weight": IntDistribution(
                min_child_low, min_child_high, log=True
            ),
        }
        return base_model, distributions

    if model_name == "svm":
        base_model = _svm_base_model(task, cuml_use_gpu)
        _, c_high, _, gamma_high = _svm_ranges(n_samples)
        distributions = {
            "model__C": FloatDistribution(0.03, c_high, log=True),
            "model__gamma": FloatDistribution(0.0003, gamma_high, log=True),
        }
        return base_model, distributions

    if model_name == "nysvm":
        ranges = _nysvm_ranges(n_samples)
        n_comp_low = ranges["n_comp_low"]
        n_comp_high = ranges["n_comp_high"]
        gamma_low = ranges["gamma_low"]
        gamma_high = ranges["gamma_high"]
        base_model = _nysvm_base_model(task, (n_comp_low + n_comp_high) // 2)
        c_key = "model__svc__C" if task == "classification" else "model__svr__C"
        distributions = {
            "model__nystroem__n_components": IntDistribution(n_comp_low, n_comp_high),
            "model__nystroem__gamma": FloatDistribution(
                gamma_low, gamma_high, log=True
            ),
            c_key: FloatDistribution(0.01, 100.0, log=True),
        }
        return base_model, distributions

    if model_name == "mlp":
        # Optuna CategoricalDistribution requires scalar choices (not tuples).
        # Encode architectures as comma-separated strings; decoded back to
        # tuples by MLPBase.fit().
        # Batch size is fixed (not tuned) and set via clone/set_params.
        hidden_presets_optuna: list[Any] = ["64", "128", "128,64", "256,128"]
        batch_sizes = _mlp_batch_sizes(n_samples)
        batch_size = batch_sizes[len(batch_sizes) // 2]
        base_model = _mlp_base_model(task)
        base_model.set_params(batch_size=batch_size)
        distributions = {
            "model__hidden_sizes": CategoricalDistribution(hidden_presets_optuna),
            "model__dropout_rate": FloatDistribution(0.20, 0.55),
            "model__learning_rate": FloatDistribution(1e-4, 3e-3, log=True),
            "model__weight_decay": FloatDistribution(1e-5, 5e-3, log=True),
        }
        return base_model, distributions

    raise ValueError(
        f"build_optuna_distributions: unsupported model '{model_name}'. "
        f"Supported: rf, xgb, svm, nysvm, mlp."
    )
