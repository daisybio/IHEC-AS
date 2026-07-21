from __future__ import annotations

"""Configuration models and constants for the splicing ML pipeline."""

import math
import os
from dataclasses import dataclass


__all__ = [
    "RNG_SEED",
    "PSI_LOW_BOUNDARY",
    "PSI_HIGH_BOUNDARY",
    "REQUIRED_COLUMNS",
    "METADATA_COLUMNS",
    "ALL_MODEL_TYPES",
    "SMOKE_MODEL_TYPES",
    "default_max_cores",
    "validate_run_config",
    "SubsetConfig",
    "RunConfig",
]


RNG_SEED = 42

# PSI boundaries shared by both tasks:
# - classification keeps extreme ranges (PSI <= low or PSI >= high)
# - regression keeps mid-range (low < PSI < high)
PSI_LOW_BOUNDARY: float = 0.2
PSI_HIGH_BOUNDARY: float = 0.8

# Required columns validated on load.
REQUIRED_COLUMNS = {
    "PSI",
    "Event Type",
    "transcript_filter",
    "Variability",
    "seqnames",
    "ontology",
}

# Non-feature columns excluded from model matrix.
METADATA_COLUMNS = {
    "ID",
    "IHEC",
    "uuid",
    "gene_id",
    "Event Type",
    "transcript_filter",
    "Variability",
    "seqnames",
    "ontology",
    "PSI",
}

# Central model registry used by parser defaults/choices and runtime validation.
ALL_MODEL_TYPES: tuple[str, ...] = (
    "linear",
    "elasticnet",
    "rf",
    "xgb",
    "lgbm",
    "beta",
    "svm",
    "nysvm",
    "mlp",
    "tabicl",
)
# Default production runs exclude:
# - svm (O(n²), unusable at SE scale)
# - nysvm (OOM-killed on SE; worst-performing model at all scales)
# - elasticnet (slower than linear, competitive only with linear/beta; not justified by performance gain)
# - rf (redundant with xgb; xgb consistently outperforms it)
# - tabicl (removed from defaults 2026-07-21, user direction) -- SE-scale runs regularly cost
#   hours-to-a-day per fold beyond xgb/lgbm/linear: its own tf_col/tf_row/tf_icl context rebuild
#   is expensive on its own, `tune_threshold_balanced_accuracy` has no tabicl exemption so it
#   refits that context 4 more times per fold for threshold tuning, and it OOM'd outright in
#   4 of 14 real logs scanned (evaluator.py's _chunked_model_call retry helps the OOMs but not
#   the underlying wall-clock cost). Still selectable explicitly via --models/include_models;
#   revisit re-adding to defaults if/once its cost profile is fixed (chunking is not the fix,
#   see the retry-loop analysis in CLAUDE.md's Known Pitfalls).
# - lgbm (added to defaults 2026-07-15, removed 2026-07-21) -- real logs across all 7
#   available classification logs confirmed xgb consistently faster (2-4x at RI scale,
#   ~13-40% at SE scale); lgbm does genuinely win on AUROC sometimes (5/5 folds outright in
#   two separate configs), so this is a real accuracy-vs-speed tradeoff, not a "never helps"
#   call -- user decided the win is too marginal to justify the added wall-clock. Still
#   selectable via --models/include_models.
# MLP was optimized with batch-size scaling, mixed precision, and early
# stopping and is now fast enough for routine production runs.
DEFAULT_MODEL_TYPES: tuple[str, ...] = tuple(
    {
        "linear",
        "xgb",
    }
)
# Smoke mode covers all available model types including MLP.
SMOKE_MODEL_TYPES: tuple[str, ...] = ALL_MODEL_TYPES


def default_max_cores() -> int:
    """Return default CPU budget (25% of available cores, min 1)."""
    total = os.cpu_count() or 1
    return max(1, math.ceil(total * 0.25))


def validate_run_config(cfg: RunConfig) -> list[str]:
    """Validate RunConfig for consistency and raise early errors.

    Returns
    -------
    list
        List of validation warnings (non-fatal).

    Raises
    ------
    ValueError
        If configuration is invalid.
    """
    warnings: list[str] = []

    # Threshold validation.
    if cfg.psi_low_threshold >= cfg.psi_high_threshold:
        raise ValueError(
            f"psi_low_threshold ({cfg.psi_low_threshold}) must be < "
            f"psi_high_threshold ({cfg.psi_high_threshold})"
        )

    if not (0.0 <= cfg.psi_low_threshold <= 1.0):
        raise ValueError(
            f"psi_low_threshold ({cfg.psi_low_threshold}) must be in [0, 1]"
        )

    if not (0.0 <= cfg.psi_high_threshold <= 1.0):
        raise ValueError(
            f"psi_high_threshold ({cfg.psi_high_threshold}) must be in [0, 1]"
        )

    # Split validation.
    if cfg.outer_splits < 2:
        raise ValueError(f"outer_splits must be >= 2, got {cfg.outer_splits}")

    if cfg.inner_splits < 2:
        raise ValueError(f"inner_splits must be >= 2, got {cfg.inner_splits}")

    if cfg.outer_splits > 20:
        warnings.append(f"outer_splits={cfg.outer_splits} is very high and may be slow")

    if cfg.only_group_col is not None and cfg.only_group_col not in {
        "seqnames",
        "ontology",
    }:
        raise ValueError(
            f"only_group_col must be one of: seqnames, ontology, got {cfg.only_group_col}"
        )

    # Task validation.
    if not cfg.run_regression and not cfg.run_classification:
        raise ValueError("Must run at least one of: regression or classification")

    # Model validation.
    allowed_models = set(ALL_MODEL_TYPES)
    invalid_models = set(cfg.include_models) - allowed_models
    if invalid_models:
        raise ValueError(f"Unknown models: {invalid_models}. Allowed: {allowed_models}")

    if len(cfg.include_models) == 0:
        raise ValueError("Must include at least one model")

    # Output validation.
    if cfg.output_level not in {"compact", "diagnostics"}:
        raise ValueError(
            f"output_level must be 'compact' or 'diagnostics', got {cfg.output_level}"
        )

    if cfg.log_level not in {"info", "debug"}:
        raise ValueError(f"log_level must be 'info' or 'debug', got {cfg.log_level}")

    # Logit regression validation.
    if not isinstance(cfg.use_logit_regression, bool):
        raise ValueError(
            f"use_logit_regression must be bool, got {type(cfg.use_logit_regression)}"
        )

    # Search strategy validation.
    if cfg.search_strategy not in {"grid", "random", "hybrid"}:
        raise ValueError(
            f"search_strategy must be one of: grid, random, hybrid, "
            f"got {cfg.search_strategy}"
        )

    # LHS scale mode validation.
    if cfg.lhs_scale_mode not in {"auto", "log", "linear"}:
        raise ValueError(
            f"lhs_scale_mode must be one of: auto, log, linear, "
            f"got {cfg.lhs_scale_mode}"
        )

    # Optuna backend validation.
    if cfg.optuna_sampler not in {"tpe", "cmaes"}:
        raise ValueError(
            f"optuna_sampler must be one of: tpe, cmaes, got {cfg.optuna_sampler}"
        )
    if cfg.optuna_n_startup_trials < 1:
        raise ValueError(
            f"optuna_n_startup_trials must be >= 1, got {cfg.optuna_n_startup_trials}"
        )
    if cfg.optuna_backend:
        try:
            import optuna  # noqa: F401
        except Exception:
            warnings.append(
                "optuna_backend=True but optuna is not installed; "
                "install with: pip install 'optuna-integration[sklearn]>=3.4'"
            )

    # Parallelism validation.
    if cfg.max_cores < 1:
        raise ValueError(f"max_cores must be >= 1, got {cfg.max_cores}")

    if cfg.param_grid_size < 1:
        raise ValueError(f"param_grid_size must be >= 1, got {cfg.param_grid_size}")

    # shap_max_samples: no lower bound to check -- <=0 is a valid sentinel
    # meaning "disable subsampling, use every test row" (publication runs).

    # Smoke mode validation.
    if cfg.smoke_mode and cfg.smoke_max_rows < 100:
        raise ValueError(f"smoke_max_rows must be >= 100, got {cfg.smoke_max_rows}")

    # Optional W&B tracking validation (non-fatal when package is missing).
    if cfg.use_wandb:
        try:
            import wandb  # noqa: F401
        except Exception:
            warnings.append(
                "use_wandb=True but wandb is not installed; falling back to NullTracker"
            )

    return warnings


@dataclass(frozen=True)
class SubsetConfig:
    """Single subset combination to evaluate."""

    event_type: str
    transcript_filter: str
    variability: str
    group_col: str


@dataclass(frozen=True)
class RunConfig:
    """Top-level runtime options controlling training and reporting."""

    data_path: str
    output_dir: str
    data_reader_backend: str = "polars"
    only_event_type: str | None = None
    only_transcript_filter: str | None = None
    only_variability: str | None = None
    only_group_col: str | None = None
    outer_splits: int = (
        5  # Original: 10; reduced to cut inner_splits 9→4 (28 vs 63 CV tasks)
    )
    inner_splits: int = 3
    random_seed: int = RNG_SEED
    max_cores: int = default_max_cores()
    param_grid_size: int = 16
    search_strategy: str = "hybrid"
    lhs_scale_mode: str = "auto"
    optuna_backend: bool = True
    optuna_sampler: str = "tpe"
    optuna_n_startup_trials: int = 5
    optuna_multivariate: bool = True
    optuna_wandb_callback: bool = True
    xgb_use_gpu: bool | None = None
    cuml_use_gpu: bool | None = None
    lgbm_use_gpu: bool | None = None
    verbose: bool = False
    include_models: tuple[str, ...] = ALL_MODEL_TYPES
    # Ablation studies: restrict model features to specific FEATURE_GROUPS
    # (see preprocessing.py), e.g. ("sequence",) for a sequence-only baseline
    # or ("histone","dnam") for epigenetics-only. None = no restriction (all
    # features, current default production behaviour).
    feature_groups: tuple[str, ...] | None = None
    run_regression: bool = True
    run_classification: bool = True
    psi_low_threshold: float = PSI_LOW_BOUNDARY
    psi_high_threshold: float = PSI_HIGH_BOUNDARY
    categorical_missing_strategy: str = "missing_token"
    categorical_missing_token: str = "__MISSING__"
    generate_html_reports: bool = True
    use_logit_regression: bool = False
    smoke_mode: bool = False
    smoke_max_rows: int = 20000
    calibrate_classifiers: bool = True
    tune_threshold: bool = True
    output_level: str = "diagnostics"
    # SHAP TreeExplainer cost scales linearly with n_test_samples (no internal
    # subsampling); on SE-scale folds (millions of rows) computing it on the
    # full outer-test fold took ~2h+ with zero progress output (found
    # 2026-07-16 via real log audit). Subsampling to a fixed size gives
    # stable mean-|SHAP| feature importance at a small, bounded cost instead.
    # <=0 disables subsampling (every test row) -- use for the final
    # publication run, where exact per-sample values matter more than
    # wall-clock. Bumped default 5000 -> 50000 (2026-07-21): real log timing
    # showed TreeExplainer at 5000 samples cost 0.4-144s even at SE scale
    # (millions of test rows) vs. xgb/lgbm's own 60-200+ *minute* searches --
    # huge unused headroom, so a bigger default gives materially more robust
    # per-feature importance for routine runs too, not just publication ones.
    shap_max_samples: int = 50000
    log_level: str = "info"
    use_wandb: bool = True
    wandb_project: str = "splicing-ml"
    wandb_entity: str | None = None
    wandb_require_auth: bool = True
    wandb_log_fold_table: bool = True
    wandb_log_tuning_details: bool = True
    wandb_log_baseline_metrics: bool = True
    wandb_fold_subruns: bool = False
