from __future__ import annotations

"""Pipeline orchestrator for splicing ML.

This module intentionally keeps only workflow orchestration logic. Most
implementation details are delegated to specialized sibling modules:
- data.py: loading/subsetting/target construction
- preprocessing.py: transformer construction
- modeling.py: tuning and evaluation
- reporting.py: per-subset HTML reports
- io_utils.py: artifact persistence
"""

import argparse
from dataclasses import asdict
from pathlib import Path
from typing import Any
import os

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline

from .config import (
    ALL_MODEL_TYPES,
    RNG_SEED,
    RunConfig,
    SMOKE_MODEL_TYPES,
    SubsetConfig,
    validate_run_config,
)
from .data import (
    benchmark_data_readers,
    build_targets,
    data_reader_startup_warnings,
    generate_subset_configs,
    load_dataset,
    subset_dataframe,
)
from .io_utils import save_json_gz, save_pickle_gz
from .metrics import bootstrap_ci, classification_metrics, regression_metrics
from .modeling import (
    balanced_group_split_indices,
    build_baseline,
    bounded_group_splits,
    evaluate_outer_fold,
    fit_best_estimator,
    metric_key,
    xgb_gpu_available,
)
from .preprocessing import add_protocol_expression_interactions, build_preprocessor
from .reporting import generate_html_reports
from .tracking import NullTracker, make_tracker
from .utils import progress_iter, set_log_level, vlog


__all__ = [
    "run_single_configuration",
    "run_all_configurations",
    "parse_args",
    "main",
]


def _configure_parallelism_guardrails(max_cores: int) -> dict[str, str]:
    """Prevent thread oversubscription by limiting BLAS/MKL threads.

    Returns
    -------
    dict
        Environment variables set (for reproducibility logging).
    """
    # Keep math libraries single-threaded so max_cores is controlled by sklearn
    # n_jobs and does not multiply with BLAS/OpenMP thread pools.
    threads_per_task = 1

    env_vars = {}
    for var in [
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        os.environ[var] = str(threads_per_task)
        env_vars[var] = str(threads_per_task)

    return env_vars


def run_single_configuration(
    df: pd.DataFrame,
    cfg: SubsetConfig,
    run_cfg: RunConfig,
    task: str,
    tracker: Any | None = None,
) -> dict[str, Any]:
    """Run nested CV for one subset+task pair.

    Steps:
    1. Subset rows by config.
    2. Build task-specific target.
    3. Build preprocessing graph.
    4. Execute grouped outer/inner CV.
    5. Aggregate fold metrics and uncertainty summary.
    """
    vlog(
        run_cfg.verbose,
        f"Starting config-task: {cfg.event_type}/{cfg.transcript_filter}/{cfg.variability}/{cfg.group_col} - {task}",
        level="info",
    )

    sdf = subset_dataframe(df, cfg, verbose=run_cfg.verbose)
    psi_response_full = sdf["PSI"].astype(float).to_numpy()
    sdf, y = build_targets(
        sdf,
        task,
        run_cfg.psi_low_threshold,
        run_cfg.psi_high_threshold,
        verbose=run_cfg.verbose,
        use_logit=(task == "regression" and run_cfg.use_logit_regression),
    )

    if run_cfg.smoke_mode and sdf.shape[0] > run_cfg.smoke_max_rows:
        pre_smoke_rows = int(sdf.shape[0])
        pre_smoke_groups = int(sdf[cfg.group_col].nunique())
        # Balanced per-outer-group sampling keeps smoke-mode representative.
        sdf = (
            sdf.groupby(cfg.group_col, group_keys=False)
            .apply(
                lambda g: g.sample(
                    n=min(
                        len(g),
                        max(
                            1,
                            run_cfg.smoke_max_rows
                            // max(1, sdf[cfg.group_col].nunique()),
                        ),
                    ),
                    random_state=run_cfg.random_seed,
                )
            )
            .reset_index(drop=True)
        )
        _, y = build_targets(
            sdf,
            task,
            run_cfg.psi_low_threshold,
            run_cfg.psi_high_threshold,
            verbose=run_cfg.verbose,
            use_logit=(task == "regression" and run_cfg.use_logit_regression),
        )
        vlog(
            run_cfg.verbose,
            "Smoke downsampling applied: "
            f"rows {pre_smoke_rows}->{int(sdf.shape[0])}, "
            f"groups={pre_smoke_groups}, target_rows={run_cfg.smoke_max_rows}",
            level="info",
        )

    # Nested-CV policy: outer_k uses one user-facing k with a minimum of 3.
    requested_outer_k = max(3, int(run_cfg.outer_splits))

    if sdf.empty or y.size < max(requested_outer_k, 10):
        return {
            "status": "skipped",
            "reason": "insufficient samples after filtering",
            "config": asdict(cfg),
            "task": task,
        }

    sdf, interaction_details = add_protocol_expression_interactions(
        sdf,
        verbose=run_cfg.verbose,
    )

    groups_outer = sdf[cfg.group_col]
    x = sdf.drop(columns=["PSI"])

    preprocessor, prep_details = build_preprocessor(
        sdf,
        categorical_missing_strategy=run_cfg.categorical_missing_strategy,
        categorical_missing_token=run_cfg.categorical_missing_token,
        verbose=run_cfg.verbose,
    )
    prep_details.update(interaction_details)

    tracker_obj = tracker if tracker is not None else NullTracker()
    tracker_obj.start_config_task(
        cfg=cfg,
        task=task,
        run_cfg=run_cfg,
        n_samples=int(sdf.shape[0]),
        prep_details=prep_details,
    )

    n_outer_splits = bounded_group_splits(
        requested_splits=requested_outer_k,
        n_groups=int(groups_outer.nunique()),
        max_splits=10,
    )
    outer_splits = balanced_group_split_indices(
        groups_outer,
        n_splits=n_outer_splits,
        seed=run_cfg.random_seed,
    )
    # Validate outer splits: ensure they partition the data without overlap.
    total_outer_idx = set()
    for tr_idx, te_idx in outer_splits:
        overlap = set(tr_idx) & set(te_idx)
        if overlap:
            raise ValueError(
                f"split_integrity_error: outer fold train/test overlap detected: "
                f"{len(overlap)} indices appear in both"
            )
        total_outer_idx.update(tr_idx)
        total_outer_idx.update(te_idx)
    if len(total_outer_idx) != len(x):
        raise ValueError(
            f"split_integrity_error: outer folds do not cover all observations "
            f"(covered {len(total_outer_idx)}/{len(x)})"
        )

    fold_results: list[dict[str, Any]] = []
    primary_scores: list[float] = []
    warnings: list[str] = []

    outer_fold_iter = progress_iter(
        list(enumerate(outer_splits, start=1)),
        total=len(outer_splits),
        desc=(
            f"folds {cfg.event_type}/{cfg.transcript_filter}/{cfg.variability}"
            f" [{task}]"
        ),
        # Keep bars in INFO verbose mode, but disable in DEBUG to avoid
        # interference with detailed per-fold log lines.
        enabled=run_cfg.verbose and run_cfg.log_level != "debug",
    )
    for fold_id, (tr_idx, te_idx) in outer_fold_iter:
        vlog(
            run_cfg.verbose,
            f"Outer fold {fold_id}: train_n={len(tr_idx)}, test_n={len(te_idx)}",
        )

        x_train = x.iloc[tr_idx]
        y_train = y[tr_idx]
        x_test = x.iloc[te_idx]
        y_test = y[te_idx]

        # Inner CV policy: use exactly the other outer folds (no re-splitting).
        local_pos = {
            int(global_i): int(local_i) for local_i, global_i in enumerate(tr_idx)
        }
        current_fold_index = fold_id - 1
        inner_splits: list[tuple[np.ndarray, np.ndarray]] = []
        inner_splits_global: list[dict[str, list[int]]] = []
        for outer_idx, (_, outer_te_idx) in enumerate(outer_splits):
            if outer_idx == current_fold_index:
                continue
            inner_valid_global = [int(i) for i in outer_te_idx if int(i) in local_pos]
            if not inner_valid_global:
                continue
            inner_valid_set = set(inner_valid_global)
            inner_valid_local = np.asarray(
                [local_pos[i] for i in inner_valid_global], dtype=int
            )
            inner_train_local = np.asarray(
                [li for li, gi in enumerate(tr_idx) if int(gi) not in inner_valid_set],
                dtype=int,
            )
            inner_splits.append((inner_train_local, inner_valid_local))
            inner_splits_global.append(
                {
                    "inner_train_indices": [
                        int(tr_idx[li]) for li in inner_train_local
                    ],
                    "inner_valid_indices": inner_valid_global,
                }
            )

        inner_splits_count = len(inner_splits)
        if inner_splits_count < 2:
            warnings.append(
                f"fold {fold_id}: insufficient remaining outer folds "
                f"for inner CV (got {inner_splits_count}); skipping fold"
            )
            continue

        vlog(
            run_cfg.verbose,
            f"Inner CV derived from outer folds: inner_splits={inner_splits_count} "
            f"(outer_k={n_outer_splits}, expected={max(0, n_outer_splits - 1)})",
        )

        # Use a single parallelism layer (inside sklearn search/calibration) to
        # keep max_cores as an effective upper bound.
        model_names = list(run_cfg.include_models)
        if task == "classification":
            model_names = [m for m in model_names if m != "beta"]

        def _fit_model_fold(model_name: str) -> dict[str, Any]:
            """Fit and evaluate one model on current fold."""
            try:
                best_estimator, tuning_info = fit_best_estimator(
                    task=task,
                    model_name=model_name,
                    x_train=x_train,
                    y_train=y_train,
                    preprocessor=preprocessor,
                    inner_cv=inner_splits,
                    max_cores=run_cfg.max_cores,
                    param_grid_size=run_cfg.param_grid_size,
                    search_strategy=run_cfg.search_strategy,
                    lhs_scale_mode=run_cfg.lhs_scale_mode,
                    xgb_use_gpu=run_cfg.xgb_use_gpu,
                    verbose=run_cfg.verbose,
                )

                ev = evaluate_outer_fold(
                    task=task,
                    best_estimator=best_estimator,
                    x_train=x_train,
                    y_train=y_train,
                    x_test=x_test,
                    y_test=y_test,
                    inner_splits=inner_splits,
                    max_cores=run_cfg.max_cores,
                    calibrate=run_cfg.calibrate_classifiers,
                    verbose=run_cfg.verbose,
                )

                # Speed optimization: skip baseline in compact mode (3-5% faster)
                if run_cfg.output_level == "compact":
                    baseline_scores = {}
                else:
                    # Baseline path mirrors the same preprocessing for fair comparison.
                    baseline_pipe = Pipeline(
                        steps=[("prep", preprocessor), ("model", build_baseline(task))]
                    )
                    baseline_pipe.fit(x_train, y_train)
                    if task == "classification":
                        baseline_prob = baseline_pipe.predict_proba(x_test)[:, 1]
                        baseline_scores = classification_metrics(
                            y_test, baseline_prob, threshold=0.5
                        )
                    else:
                        baseline_pred = baseline_pipe.predict(x_test)
                        baseline_scores = regression_metrics(y_test, baseline_pred)

                pkey = metric_key(task)
                pscore = float(ev["scores"][pkey])

                fold_result = {
                    "model_name": model_name,
                    "scores": ev["scores"],
                    "threshold": ev["threshold"],
                    "tuning": tuning_info,
                    "baseline_scores": baseline_scores,
                    "delta_vs_baseline": {
                        k: float(
                            ev["scores"].get(k, np.nan) - baseline_scores.get(k, np.nan)
                        )
                        for k in ev["scores"].keys()
                    },
                    "outer_fold": fold_id,
                    "outer_test_indices": te_idx.tolist(),
                    "outer_test_groups": groups_outer.iloc[te_idx].tolist(),
                    "inner_splits_global": inner_splits_global,
                    "primary_score": pscore,
                }
                # Include predictions only in diagnostics mode.
                if run_cfg.output_level == "diagnostics":
                    fold_result["y_true"] = ev["y_true"].tolist()
                    fold_result["y_pred"] = ev["y_pred"].tolist()
                fold_result["_estimator"] = ev.get("estimator")

                return fold_result
            except Exception as exc:
                return {"error": str(exc), "model_name": model_name, "fold_id": fold_id}

        model_results = [_fit_model_fold(m) for m in model_names]

        # Process results and collect warnings
        for result in model_results:
            if "error" in result:
                warn_msg = (
                    f"fold {result['fold_id']} model {result['model_name']}: "
                    f"{result['error']}"
                )
                warnings.append(warn_msg)
                vlog(run_cfg.verbose, f"Model failed: {warn_msg}", level="info")
            else:
                fold_results.append(result)
                estimator = result.pop("_estimator", None)
                pscore = result.pop("primary_score")
                primary_scores.append(pscore)
                tracker_obj.log_fold_result(
                    fold_id=int(fold_id),
                    model_name=str(result["model_name"]),
                    scores=result.get("scores", {}),
                    tuning_info=result.get("tuning", {}),
                    y_true=result.get("y_true", []),
                    y_pred=result.get("y_pred", []),
                    estimator=estimator,
                    task=task,
                    x_train=x_train,
                    x_test=x_test,
                )
                vlog(
                    run_cfg.verbose,
                    f"Completed fold={fold_id}, model={result['model_name']}, primary_metric={metric_key(task)}={pscore:.6f}",
                )

    if not fold_results:
        fail_result = {
            "status": "failed",
            "reason": "no successful fold/model runs",
            "warnings": warnings,
            "config": asdict(cfg),
            "task": task,
        }
        tracker_obj.log_config_task_summary(fail_result)
        tracker_obj.finish_config_task()
        return {
            "status": "failed",
            "reason": "no successful fold/model runs",
            "warnings": warnings,
            "config": asdict(cfg),
            "task": task,
        }

    ci_low, ci_high = bootstrap_ci(np.asarray(primary_scores), seed=run_cfg.random_seed)
    out = {
        "status": "ok",
        "config": asdict(cfg),
        "task": task,
        "n_samples": int(sdf.shape[0]),
        "response_distribution": {
            # Keep compact but informative payload for plotting/reporting.
            "psi_sample": psi_response_full[:20000].tolist(),
            "psi_sample_size": int(min(20000, psi_response_full.shape[0])),
            "psi_total_size": int(psi_response_full.shape[0]),
            "binarization_thresholds": (
                [
                    float(run_cfg.psi_low_threshold),
                    float(run_cfg.psi_high_threshold),
                ]
                if task == "classification"
                else None
            ),
            "label_distribution": (
                {
                    "0": int(np.sum(y == 0)) if task == "classification" else None,
                    "1": int(np.sum(y == 1)) if task == "classification" else None,
                }
                if task == "classification"
                else None
            ),
        },
        "preprocessing": prep_details,
        "grouping": {
            "group_col": cfg.group_col,
        },
        "primary_metric": metric_key(task),
        "primary_metric_mean": float(np.mean(primary_scores)),
        "primary_metric_std": float(np.std(primary_scores, ddof=0)),
        "primary_metric_bootstrap_ci": [ci_low, ci_high],
        "fold_results": fold_results,
        "warnings": warnings,
        "reproducibility": {
            "random_seed": run_cfg.random_seed,
            "psi_thresholds": [run_cfg.psi_low_threshold, run_cfg.psi_high_threshold],
        },
    }
    vlog(
        run_cfg.verbose,
        f"Finished config-task with status=ok, n_fold_results={len(fold_results)}",
        level="info",
    )
    tracker_obj.log_config_task_summary(out)
    tracker_obj.finish_config_task()
    return out


def run_all_configurations(run_cfg: RunConfig) -> dict[str, Any]:
    """Run all subset/task combinations and persist outputs."""
    np.random.seed(run_cfg.random_seed)
    vlog(run_cfg.verbose, "Starting full run_all_configurations", level="info")
    tracker = make_tracker(run_cfg)
    tracker.start_run(run_cfg)

    # Speed optimization: prefer polars for loading (20-40% faster) unless explicitly disabled
    backend = run_cfg.data_reader_backend
    if backend == "auto":
        try:
            import polars as pl  # noqa: F401

            backend = "polars"
        except Exception:
            backend = "pandas"

    df = load_dataset(
        run_cfg.data_path,
        verbose=run_cfg.verbose,
        reader_backend=backend,
    )
    configs = generate_subset_configs(df, verbose=run_cfg.verbose)
    if run_cfg.only_event_type is not None:
        configs = [c for c in configs if c.event_type == run_cfg.only_event_type]
    if run_cfg.only_transcript_filter is not None:
        configs = [
            c for c in configs if c.transcript_filter == run_cfg.only_transcript_filter
        ]
    if run_cfg.only_variability is not None:
        configs = [c for c in configs if c.variability == run_cfg.only_variability]
    if run_cfg.only_group_col is not None:
        configs = [c for c in configs if c.group_col == run_cfg.only_group_col]
    vlog(run_cfg.verbose, f"Configs after CLI filters: {len(configs)}", level="info")
    if run_cfg.smoke_mode:
        configs = configs[:1]
        vlog(
            run_cfg.verbose,
            "Smoke mode enabled: limiting to first subset config",
            level="info",
        )

    out_dir = Path(run_cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks: list[str] = []
    if run_cfg.run_regression:
        tasks.append("regression")
    if run_cfg.run_classification:
        tasks.append("classification")

    results: list[dict[str, Any]] = []
    work_items = [(cfg, task) for cfg in configs for task in tasks]
    for cfg, task in work_items:
        # Add a visual separator so tqdm task bars do not appear merged.
        vlog(
            run_cfg.verbose,
            f"--- Starting {task} task for config {cfg} ---",
            level="info",
        )
        results.append(
            run_single_configuration(df, cfg, run_cfg, task, tracker=tracker)
        )

    if run_cfg.generate_html_reports:
        generate_html_reports(results, out_dir, verbose=run_cfg.verbose)

    bundle = {"run_config": asdict(run_cfg), "results": results}
    pkl_path = out_dir / "splicing_ml_results.pkl.gz"
    json_path = out_dir / "splicing_ml_results.json.gz"
    save_pickle_gz(bundle, pkl_path)
    save_json_gz(bundle, json_path)
    tracker.log_artifact(
        pkl_path, name="splicing-ml-results-pkl", artifact_type="results"
    )
    tracker.log_artifact(
        json_path, name="splicing-ml-results-json", artifact_type="results"
    )
    tracker.finish_run(results)
    vlog(run_cfg.verbose, f"Artifacts written to {out_dir}", level="info")
    return bundle


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for end-to-end execution."""
    from .config import default_max_cores

    p = argparse.ArgumentParser(description="Nested CV splicing ML pipeline")
    p.add_argument(
        "--data-path", default="processed_data/aggregated_dt_filtered.csv.gz"
    )
    p.add_argument(
        "--data-reader-backend",
        default="auto",
        choices=["auto", "pandas", "polars"],
    )
    p.add_argument(
        "--benchmark-data-readers",
        action="store_true",
        help="Benchmark auto/pandas/polars loading backends on --data-path and exit",
    )
    p.add_argument(
        "--benchmark-runs",
        type=int,
        default=1,
        help="Number of repeated runs per backend for --benchmark-data-readers",
    )
    p.add_argument("--output-dir", default="splicing_ml/output/ml_splicing_outputs")
    p.add_argument("--only-event-type", choices=["SE", "RI"], default=None)
    p.add_argument(
        "--only-transcript-filter",
        choices=["transcripts", "tsl_filtered", "biotype_filtered"],
        default=None,
    )
    p.add_argument(
        "--only-variability",
        choices=["High", "Low", "both"],
        default=None,
    )
    p.add_argument(
        "--only-group",
        default=None,
        choices=["seqnames", "ontology"],
        help="Optional filter to run only one grouping variable; default runs both.",
    )
    p.add_argument(
        "--models",
        nargs="+",
        default=list(ALL_MODEL_TYPES),
        choices=list(ALL_MODEL_TYPES),
    )
    p.add_argument("--seed", type=int, default=RNG_SEED)
    p.add_argument("--max-cores", type=int, default=default_max_cores())
    p.add_argument("--param-grid-size", type=int, default=16)
    p.add_argument(
        "--search-strategy",
        type=str,
        default="hybrid",
        choices=["grid", "random", "hybrid"],
    )
    p.add_argument(
        "--lhs-scale-mode",
        type=str,
        default="auto",
        choices=["auto", "log", "linear"],
        help="Scale mapping for LHS candidates when search uses random/hybrid",
    )
    p.add_argument("--outer-splits", type=int, default=10)
    p.add_argument("--skip-regression", action="store_true")
    p.add_argument("--skip-classification", action="store_true")
    p.add_argument("--verbose", action="store_true")
    p.add_argument(
        "--debug",
        action="store_true",
        help="Show detailed debug logs (implies --verbose)",
    )
    p.add_argument(
        "--categorical-missing-strategy",
        default="missing_token",
        choices=["missing_token", "most_frequent"],
    )
    p.add_argument("--categorical-missing-token", default="__MISSING__")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--smoke-max-rows", type=int, default=20000)
    p.add_argument("--no-html-reports", action="store_true")
    p.add_argument(
        "--calibration",
        action="store_true",
        help="Enable classifier calibration (disabled by default for speed)",
    )
    p.add_argument(
        "--output-level",
        type=str,
        default="diagnostics",
        choices=["compact", "diagnostics"],
        help="'compact' stores aggregated metrics only; 'diagnostics' includes per-fold predictions",
    )
    p.add_argument(
        "--disable-logit-regression",
        action="store_true",
        help="Disable logit transform for regression targets (model PSI on raw [0,1] scale)",
    )
    p.add_argument("--wandb", action="store_true")
    p.add_argument("--wandb-project", default="splicing-ml")
    p.add_argument("--wandb-entity", default=None)
    return p.parse_args()


def main() -> None:
    """CLI entrypoint."""
    args = parse_args()

    if args.debug:
        args.verbose = True

    set_log_level("debug" if args.debug else ("info" if args.verbose else "none"))

    if args.smoke:
        args.data_path = "processed_data/aggregated_dt_filtered.validation300k.csv.gz"
        # Keep smoke runs lightweight regardless of user default grid size.
        args.param_grid_size = min(args.param_grid_size, 8)

    startup_warnings = data_reader_startup_warnings(args.data_reader_backend)
    for w in startup_warnings:
        vlog(args.verbose, f"Warning: {w}", level="info")

    if args.benchmark_data_readers:
        bench = benchmark_data_readers(
            args.data_path,
            runs=max(1, args.benchmark_runs),
            verbose=args.verbose,
        )
        vlog(args.verbose, "Data reader benchmark summary", level="info")
        vlog(args.verbose, f"Path: {bench['path']}", level="info")
        vlog(args.verbose, f"Runs per backend: {bench['runs']}", level="info")
        caps = bench["capabilities"]
        vlog(
            args.verbose,
            f"Capabilities: polars={caps['polars']}, pyarrow={caps['pyarrow']}",
            level="info",
        )
        for row in bench["results"]:
            if row["ok"]:
                vlog(
                    args.verbose,
                    " - "
                    f"{row['backend']}: mean={row['mean_seconds']:.3f}s "
                    f"(min={row['min_seconds']:.3f}s, max={row['max_seconds']:.3f}s), "
                    f"shape={tuple(row['shape'])}",
                    level="info",
                )
            else:
                vlog(
                    args.verbose,
                    f" - {row['backend']}: FAILED ({row['error']})",
                    level="info",
                )
        if bench["winner_backend"]:
            vlog(
                args.verbose,
                f"Recommended backend: {bench['winner_backend']}",
                level="info",
            )
        else:
            vlog(args.verbose, "No backend succeeded", level="info")
        return

    models = tuple(args.models)
    if args.smoke:
        # extract all available models in smoke mode for comprehensive testing, regardless of user default
        models = SMOKE_MODEL_TYPES

    xgb_use_gpu = None
    if "xgb" in models:
        xgb_use_gpu = xgb_gpu_available()
        vlog(
            args.verbose,
            f"XGBoost runtime device check: {'cuda' if xgb_use_gpu else 'cpu'}",
            level="info",
        )

    outer_splits = args.outer_splits
    if args.smoke:
        outer_splits = 3  # Exactly 3 outer folds in smoke mode for robust evaluation
    # Inner CV is always derived from outer folds: one inner fold per remaining outer fold.
    inner_splits = max(2, outer_splits - 1)

    run_cfg = RunConfig(
        data_path=args.data_path,
        output_dir=args.output_dir,
        data_reader_backend=args.data_reader_backend,
        only_event_type=args.only_event_type,
        only_transcript_filter=args.only_transcript_filter,
        only_variability=args.only_variability,
        only_group_col=args.only_group,
        outer_splits=outer_splits,
        inner_splits=inner_splits,
        random_seed=args.seed,
        max_cores=max(1, args.max_cores),
        param_grid_size=max(1, args.param_grid_size),
        search_strategy=args.search_strategy,
        lhs_scale_mode=args.lhs_scale_mode,
        xgb_use_gpu=xgb_use_gpu,
        verbose=args.verbose,
        include_models=models,
        run_regression=not args.skip_regression,
        run_classification=not args.skip_classification,
        use_logit_regression=not args.disable_logit_regression,
        categorical_missing_strategy=args.categorical_missing_strategy,
        categorical_missing_token=args.categorical_missing_token,
        generate_html_reports=not args.no_html_reports,
        smoke_mode=args.smoke,
        smoke_max_rows=args.smoke_max_rows,
        calibrate_classifiers=args.calibration,
        output_level=args.output_level,
        log_level=("debug" if args.debug else "info"),
        use_wandb=args.wandb,
        wandb_project=args.wandb_project,
        wandb_entity=args.wandb_entity,
    )

    # Validate configuration early.
    val_warnings = validate_run_config(run_cfg)
    for w in val_warnings:
        vlog(args.verbose, f"Config warning: {w}", level="info")

    vlog(
        args.verbose,
        f"Running with max_cores={run_cfg.max_cores}, "
        f"param_grid_size={run_cfg.param_grid_size}, "
        f"search_strategy={run_cfg.search_strategy}, "
        f"lhs_scale_mode={run_cfg.lhs_scale_mode}, "
        f"xgb_device={'cuda' if run_cfg.xgb_use_gpu else 'cpu'}, "
        f"data_reader_backend={run_cfg.data_reader_backend}, "
        f"data_path={run_cfg.data_path}, "
        f"group_filter={run_cfg.only_group_col}, "
        f"logit_regression={run_cfg.use_logit_regression}, "
        f"outer_splits={run_cfg.outer_splits}, inner_splits={run_cfg.inner_splits}, "
        f"calibrate_classifiers={run_cfg.calibrate_classifiers}, "
        f"output_level={run_cfg.output_level}, "
        f"log_level={run_cfg.log_level}",
        level="info",
    )

    # Configure parallelism guardrails to prevent thread oversubscription.
    thread_env = _configure_parallelism_guardrails(run_cfg.max_cores)
    vlog(args.verbose, f"Parallelism guardrails: {thread_env}", level="info")

    bundle = run_all_configurations(run_cfg)
    vlog(
        args.verbose,
        f"Completed {len(bundle['results'])} config-task runs",
        level="info",
    )


if __name__ == "__main__":
    main()
