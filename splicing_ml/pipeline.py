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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import os

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline

from .config import (
    ALL_MODEL_TYPES,
    DEFAULT_MODEL_TYPES,
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
    estimate_filter_fraction_for_path,
    generate_subset_configs,
    load_dataset,
    recommend_auto_data_reader_backend,
    recommend_polars_max_threads,
    subset_dataframe,
)
from .io_utils import save_json_gz, save_pickle_gz
from .metrics import bootstrap_ci, classification_metrics, regression_metrics
from .modeling import (
    balanced_group_split_indices,
    build_baseline,
    bounded_group_splits,
    cuml_gpu_available,
    evaluate_outer_fold,
    fit_best_estimator,
    metric_key,
    xgb_gpu_available,
)
from .ontology_hierarchy import map_ontology_to_supergroups
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

    DEPRECATED: when running under SLURM, cgroup limits (set by -c N) already
    constrain available CPUs, so this explicit env-var approach is redundant.
    Retained for reference in case of non-SLURM deployments.

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

    # Nested-CV policy: outer_k uses one user-facing k with a minimum of 4.
    requested_outer_k = max(4, int(run_cfg.outer_splits))

    grouping_col = cfg.group_col
    grouping_details: dict[str, Any] = {
        "group_col": cfg.group_col,
        "grouping_col_resolved": cfg.group_col,
        "supergrouping_enabled": False,
    }
    if cfg.group_col == "ontology":
        unique_ontology_n = int(
            sdf["ontology"].astype(str).str.strip().str.lower().nunique()
        )
        ontology_requested_outer_k = max(4, requested_outer_k)
        if unique_ontology_n >= 4:
            target_supergroups = min(ontology_requested_outer_k, unique_ontology_n)
            ontology_supergroups, hierarchy_diag = map_ontology_to_supergroups(
                sdf["ontology"],
                n_groups=target_supergroups,
            )
            sdf = sdf.assign(__group_outer=ontology_supergroups.to_numpy())
            grouping_col = "__group_outer"
            grouping_details.update(hierarchy_diag)
            grouping_details["grouping_col_resolved"] = grouping_col
            vlog(
                run_cfg.verbose,
                "Using dynamic hierarchical ontology grouping "
                f"(requested_outer_splits={requested_outer_k}, ontology_requested_outer_splits={ontology_requested_outer_k}, unique_ontology={unique_ontology_n}, "
                f"target_supergroups={target_supergroups})",
                level="info",
            )
        else:
            grouping_details["supergrouping_scheme"] = "original"
            grouping_details["supergrouping_scheme_kind"] = (
                "insufficient_unique_ontology"
            )
            grouping_details["original_ontology_total"] = unique_ontology_n
            vlog(
                run_cfg.verbose,
                "Using original ontology groups: fewer than 4 unique ontology labels are available "
                f"(unique_ontology={unique_ontology_n})",
                level="info",
            )

    if run_cfg.smoke_mode and sdf.shape[0] > run_cfg.smoke_max_rows:
        pre_smoke_rows = int(sdf.shape[0])
        pre_smoke_groups = int(sdf[grouping_col].nunique())
        # Balanced per-outer-group sampling keeps smoke-mode representative.
        per_group_target = max(
            1,
            run_cfg.smoke_max_rows // max(1, pre_smoke_groups),
        )
        sampled_parts: list[pd.DataFrame] = []
        for group_name, group_df in sdf.groupby(grouping_col, sort=False):
            group_df_local = group_df
            if grouping_col not in group_df_local.columns:
                # pandas groupby behavior changed across versions; keep the key.
                group_df_local = group_df_local.assign(**{grouping_col: group_name})
            sampled_parts.append(
                group_df_local.sample(
                    n=min(len(group_df_local), per_group_target),
                    random_state=run_cfg.random_seed,
                )
            )
        sdf = pd.concat(sampled_parts, ignore_index=True)
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

    groups_outer = sdf[grouping_col]
    feature_drop_columns = ["PSI"]
    if grouping_col not in {"seqnames", "ontology"} and grouping_col in sdf.columns:
        # Keep derived grouping helpers out of model features.
        feature_drop_columns.append(grouping_col)
    x = sdf.drop(columns=feature_drop_columns)

    preprocessor, prep_details = build_preprocessor(
        x,
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
        y=y if task == "classification" else None,
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
        outer_test_group_set = set(groups_outer.iloc[te_idx].astype(str).tolist())

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

            # Split-integrity checks: validation groups must be disjoint from
            # both inner-train groups and current outer-test groups.
            inner_train_global = tr_idx[inner_train_local]
            inner_valid_global_arr = tr_idx[inner_valid_local]
            inner_train_groups = set(
                groups_outer.iloc[inner_train_global].astype(str).tolist()
            )
            inner_valid_groups = set(
                groups_outer.iloc[inner_valid_global_arr].astype(str).tolist()
            )
            tv_overlap = inner_train_groups & inner_valid_groups
            if tv_overlap:
                raise ValueError(
                    "split_integrity_error: inner train/validation group overlap "
                    f"in outer fold {fold_id} (overlap_groups={sorted(tv_overlap)})"
                )
            test_val_overlap = outer_test_group_set & inner_valid_groups
            if test_val_overlap:
                raise ValueError(
                    "split_integrity_error: inner validation/outer test group overlap "
                    f"in outer fold {fold_id} (overlap_groups={sorted(test_val_overlap)})"
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
        if run_cfg.log_level == "debug":

            def _grp(idx_arr) -> str:
                """Group-name:count summary for a list/array of global row indices."""
                s = groups_outer.iloc[idx_arr].astype(str).value_counts().sort_index()
                return ", ".join(f"{g}:{int(n)}" for g, n in s.items()) or "<empty>"

            vlog(
                run_cfg.verbose,
                f"Fold {fold_id} outer TEST  (n={len(te_idx)}): {_grp(te_idx)}",
                level="debug",
            )
            vlog(
                run_cfg.verbose,
                f"Fold {fold_id} outer TRAIN (n={len(tr_idx)}): {_grp(tr_idx)}",
                level="debug",
            )
            n_inner = len(inner_splits_global)
            # For each split, ES val = smallest scoring-val from any other split.
            # ES train = split's train minus that borrowed val.
            es_val_source: list[int] = []  # 1-based split number of the borrowed val
            for i in range(n_inner):
                es_src = min(
                    (j for j in range(n_inner) if j != i),
                    key=lambda j: len(inner_splits_global[j]["inner_valid_indices"]),
                )
                es_val_source.append(es_src + 1)  # convert to 1-based

            for split_num, split_info in enumerate(inner_splits_global, start=1):
                scoring_val_idx = split_info["inner_valid_indices"]
                full_train_idx = split_info["inner_train_indices"]
                # Borrow the smallest other split's val as ES hold-out
                es_src_num = es_val_source[split_num - 1]
                es_val_idx = inner_splits_global[es_src_num - 1]["inner_valid_indices"]
                es_val_set = set(es_val_idx)
                es_train_idx = [i for i in full_train_idx if i not in es_val_set]
                vlog(
                    run_cfg.verbose,
                    f"Fold {fold_id} inner split{split_num} "
                    f"SCORING_VAL (n={len(scoring_val_idx)}): {_grp(scoring_val_idx)}",
                    level="debug",
                )
                vlog(
                    run_cfg.verbose,
                    f"Fold {fold_id} inner split{split_num} "
                    f"ES_VAL      (n={len(es_val_idx)}, from split{es_src_num}): "
                    f"{_grp(es_val_idx)}",
                    level="debug",
                )
                vlog(
                    run_cfg.verbose,
                    f"Fold {fold_id} inner split{split_num} "
                    f"ES_TRAIN    (n={len(es_train_idx)}): {_grp(es_train_idx)}",
                    level="debug",
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
                    cuml_use_gpu=run_cfg.cuml_use_gpu,
                    verbose=run_cfg.verbose,
                    debug_grid_progress=run_cfg.log_level == "debug",
                    optuna_backend=run_cfg.optuna_backend,
                    optuna_sampler=run_cfg.optuna_sampler,
                    optuna_n_startup_trials=run_cfg.optuna_n_startup_trials,
                    optuna_multivariate=run_cfg.optuna_multivariate,
                    optuna_wandb_callback=run_cfg.optuna_wandb_callback,
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
                    tune_threshold=run_cfg.tune_threshold,
                    verbose=run_cfg.verbose,
                    psi_low=run_cfg.psi_low_threshold,
                    psi_high=run_cfg.psi_high_threshold,
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
                        # Mirror evaluate_outer_fold dual-scale logic so baseline_scores
                        # has original_/logit_ prefixed keys matching regular model scores.
                        from .models.beta import _inverse_logit, _logit_transform

                        _bp = np.asarray(baseline_pred, dtype=float)
                        _yt = np.asarray(y_test, dtype=float)
                        _is_logit = bool(np.any(_yt < 0.0) or np.any(_yt > 1.0))
                        if _is_logit:
                            _yt_orig = _inverse_logit(_yt)
                            _bp_orig = _inverse_logit(_bp)
                        else:
                            _yt_orig = _yt
                            _bp_orig = _bp
                        _scores_orig = regression_metrics(_yt_orig, _bp_orig)
                        _yt_logit = _logit_transform(_yt_orig)
                        _bp_logit = _logit_transform(_bp_orig)
                        _scores_logit = regression_metrics(_yt_logit, _bp_logit)
                        baseline_scores = dict(_scores_orig)
                        for _k, _v in _scores_orig.items():
                            baseline_scores[f"original_{_k}"] = float(_v)
                        for _k, _v in _scores_logit.items():
                            baseline_scores[f"logit_{_k}"] = float(_v)

                # Compute training scores to detect overfitting.
                fitted_estimator = ev.get("estimator")
                train_scores: dict[str, float] = {}
                if fitted_estimator is not None and run_cfg.output_level == "diagnostics":
                    try:
                        if task == "classification":
                            if hasattr(fitted_estimator, "predict_proba"):
                                y_train_prob = fitted_estimator.predict_proba(x_train)[
                                    :, 1
                                ]
                            else:
                                _decision = fitted_estimator.decision_function(x_train)
                                y_train_prob = 1.0 / (1.0 + np.exp(-_decision))

                            train_scores = classification_metrics(
                                y_train,
                                y_train_prob,
                                threshold=ev.get("threshold", 0.5),
                            )
                        else:
                            from .models.beta import _inverse_logit, _logit_transform

                            y_train_arr = np.asarray(y_train, dtype=float)
                            y_train_pred = np.asarray(
                                fitted_estimator.predict(x_train), dtype=float
                            )
                            # Mirror evaluator scale logic: convert to PSI scale if needed.
                            is_logit = bool(
                                np.any(y_train_arr < 0.0) or np.any(y_train_arr > 1.0)
                            )
                            if is_logit:
                                y_train_orig = _inverse_logit(y_train_arr)
                                from .models.evaluator import _model_outputs_psi_scale

                                y_train_pred_orig = (
                                    np.clip(y_train_pred, 0.0, 1.0)
                                    if _model_outputs_psi_scale(fitted_estimator)
                                    else _inverse_logit(y_train_pred)
                                )
                            else:
                                y_train_orig = y_train_arr
                                y_train_pred_orig = y_train_pred
                            train_scores = regression_metrics(
                                y_train_orig, y_train_pred_orig
                            )
                    except Exception:
                        pass

                pkey = metric_key(task)
                pscore = float(ev["scores"][pkey])

                fold_result = {
                    "model_name": model_name,
                    "scores": ev["scores"],
                    "train_scores": train_scores,
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
                    "psi_bin_metrics": ev.get("psi_bin_metrics", {}),
                }
                # Include predictions only in diagnostics mode.
                if run_cfg.output_level == "diagnostics":
                    fold_result["y_true"] = ev["y_true"].tolist()
                    fold_result["y_pred"] = ev["y_pred"].tolist()
                    if task == "classification":
                        fold_result["baseline_y_pred"] = baseline_prob.tolist()
                    else:
                        fold_result["baseline_y_pred"] = baseline_pred.tolist()
                fold_result["_estimator"] = ev.get("estimator")

                # Release TabICL's training-data cache and GPU/CPU tensors
                # immediately so subsequent folds don't accumulate RAM.
                if model_name == "tabicl":
                    del best_estimator, fitted_estimator
                    try:
                        import torch as _t
                        if _t.cuda.is_available():
                            _t.cuda.empty_cache()
                    except Exception:
                        pass

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
                    train_scores=result.get("train_scores", {}),
                    primary_metric=metric_key(task),
                    baseline_scores=result.get("baseline_scores", {}),
                    delta_vs_baseline=result.get("delta_vs_baseline", {}),
                )
                train_primary = result.get("train_scores", {}).get(metric_key(task))
                train_str = (
                    f", train_{metric_key(task)}={train_primary:.6f}"
                    if train_primary is not None
                    else ""
                )
                vlog(
                    run_cfg.verbose,
                    f"Completed fold={fold_id}, model={result['model_name']}, primary_metric={metric_key(task)}={pscore:.6f}{train_str}",
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
            "binarization_thresholds": [
                float(run_cfg.psi_low_threshold),
                float(run_cfg.psi_high_threshold),
            ],
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
            **grouping_details,
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
    run_datetime = datetime.now(timezone.utc).isoformat()
    np.random.seed(run_cfg.random_seed)
    vlog(run_cfg.verbose, "Starting full run_all_configurations", level="info")
    tracker = make_tracker(run_cfg)
    tracker.start_run(run_cfg)

    # Choose loader backend with memory-aware guardrails.
    backend = run_cfg.data_reader_backend
    if backend == "auto":
        backend, backend_reason = recommend_auto_data_reader_backend(
            path=run_cfg.data_path,
            max_cores=run_cfg.max_cores,
            filter_event_type=run_cfg.only_event_type,
            filter_transcript_filter=run_cfg.only_transcript_filter,
            filter_variability=run_cfg.only_variability,
        )
        vlog(
            run_cfg.verbose,
            f"Auto data reader selected {backend}: {backend_reason}",
            level="info",
        )

    df = load_dataset(
        run_cfg.data_path,
        verbose=run_cfg.verbose,
        reader_backend=backend,
        filter_event_type=run_cfg.only_event_type,
        filter_transcript_filter=run_cfg.only_transcript_filter,
        filter_variability=run_cfg.only_variability,
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
        generate_html_reports(
            results, out_dir, verbose=run_cfg.verbose, run_datetime=run_datetime
        )

    # Save per-task artifacts so classification and regression can be loaded
    # and reported independently.
    results_by_task: dict[str, list[dict[str, Any]]] = {
        "classification": [],
        "regression": [],
    }
    for r in results:
        t = r.get("task", "")
        if t in results_by_task:
            results_by_task[t].append(r)

    for task_name, task_results in results_by_task.items():
        if not task_results:
            continue
        task_bundle = {
            "run_datetime": run_datetime,
            "run_config": asdict(run_cfg),
            "results": task_results,
        }
        pkl_path = out_dir / f"splicing_ml_results_{task_name}.pkl.gz"
        json_path = out_dir / f"splicing_ml_results_{task_name}.json.gz"
        save_pickle_gz(task_bundle, pkl_path)
        save_json_gz(task_bundle, json_path)
        tracker.log_artifact(
            pkl_path,
            name=f"splicing-ml-results-{task_name}-pkl",
            artifact_type="results",
        )
        tracker.log_artifact(
            json_path,
            name=f"splicing-ml-results-{task_name}-json",
            artifact_type="results",
        )
        vlog(
            run_cfg.verbose,
            f"Saved {task_name} artifacts: {pkl_path.name}, {json_path.name}",
            level="info",
        )

    bundle = {
        "run_datetime": run_datetime,
        "run_config": asdict(run_cfg),
        "results": results,
    }
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
        default="polars",
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
        default="seqnames",
        choices=["seqnames", "ontology"],
        help="Filter on which grouping variable to use; default runs seqnames/chromosomes (unseen events).",
    )
    p.add_argument(
        "--models",
        nargs="+",
        default=list(DEFAULT_MODEL_TYPES),
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
    p.add_argument(
        "--optuna",
        action="store_true",
        help="Use Optuna Bayesian optimisation (TPE/CMA-ES) instead of GridSearchCV",
    )
    p.add_argument(
        "--optuna-sampler",
        type=str,
        default="tpe",
        choices=["tpe", "cmaes"],
        help="Optuna sampler: 'tpe' (Tree-structured Parzen Estimator) or 'cmaes'",
    )
    p.add_argument(
        "--optuna-n-startup-trials",
        type=int,
        default=5,
        help="Number of random trials before Optuna's surrogate model takes over",
    )
    p.add_argument(
        "--optuna-wandb-callback",
        action="store_true",
        help="Log per-trial Optuna metrics to the active W&B run (requires --wandb)",
    )
    # Reduced from 10 → 5: inner_splits = outer_splits - 1, so 10 outer folds
    # gives 9 inner folds → 63 LogisticRegressionCV tasks vs 28 at 5 outer folds.
    # Original: default=10
    p.add_argument("--outer-splits", type=int, default=5)
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
        "--no-tune-threshold",
        dest="tune_threshold",
        action="store_false",
        help="Disable decision threshold tuning on inner-fold predictions (default: enabled)",
    )
    p.set_defaults(tune_threshold=True)
    p.add_argument(
        "--output-level",
        type=str,
        default="diagnostics",
        choices=["compact", "diagnostics"],
        help="'compact' stores aggregated metrics only; 'diagnostics' includes per-fold predictions",
    )
    p.add_argument(
        "--logit-regression",
        action="store_true",
        help="Apply logit transform to regression targets (model on log-odds scale instead of raw PSI)",
    )
    p.add_argument("--wandb", action="store_true")
    p.add_argument("--wandb-project", default="splicing-ml")
    p.add_argument("--wandb-entity", default=None)
    p.add_argument(
        "--wandb-no-require-auth",
        action="store_true",
        help="Allow fallback to NullTracker when W&B auth is missing or invalid",
    )
    p.add_argument(
        "--wandb-no-fold-table",
        action="store_true",
        help="Disable fold-level comparison tables in W&B child runs",
    )
    p.add_argument(
        "--wandb-no-tuning-details",
        action="store_true",
        help="Disable structured tuning detail logging in W&B",
    )
    p.add_argument(
        "--wandb-no-baseline-metrics",
        action="store_true",
        help="Disable baseline and delta metric logging in W&B",
    )
    p.add_argument(
        "--wandb-fold-subruns",
        action="store_true",
        help="Create optional per-fold W&B sub-runs (higher run/API volume)",
    )
    p.add_argument(
        "--check-gpu",
        action="store_true",
        help="Print GPU availability for each backend and exit",
    )
    return p.parse_args()


def _run_gpu_check() -> None:
    """Print GPU availability for XGBoost, cuML, and PyTorch, then exit."""
    # XGBoost
    try:
        from splicing_ml.models.xgb_utils import xgb_gpu_available

        xgb_ok = xgb_gpu_available()
        print(f"XGBoost : {'cuda' if xgb_ok else 'cpu (no GPU)'}")
    except Exception as e:
        print(f"XGBoost : error ({e})")

    # cuML SVM
    try:
        from splicing_ml.models.cuml_utils import cuml_gpu_available

        cuml_ok = cuml_gpu_available()
        print(f"cuML SVM: {'cuda' if cuml_ok else 'cpu (no GPU)'}")
    except Exception as e:
        print(f"cuML SVM: error ({e})")

    # PyTorch / MLP
    try:
        import torch as _torch

        device = "cpu"
        try:
            if _torch.accelerator.is_available():
                device = str(_torch.accelerator.current_accelerator())
        except AttributeError:
            if _torch.cuda.is_available():
                device = "cuda"
            elif hasattr(_torch.backends, "mps") and _torch.backends.mps.is_available():
                device = "mps"
        print(f"PyTorch : {device}")
    except ImportError:
        print("PyTorch : not installed")


def main() -> None:
    """CLI entrypoint."""
    args = parse_args()

    if args.debug:
        args.verbose = True

    set_log_level("debug" if args.debug else ("info" if args.verbose else "none"))

    if args.check_gpu:
        _run_gpu_check()
        return

    if args.smoke:
        args.data_path = "processed_data/aggregated_dt_filtered.validation300k.csv.gz"
        # Keep smoke runs lightweight regardless of user default grid size.
        args.param_grid_size = min(args.param_grid_size, 4)
        # Always exercise the Optuna path in smoke so both search backends are covered.
        args.optuna = True

    requested_models = tuple(args.models)

    vlog(
        args.verbose,
        f"Requested run args: max_cores={max(1, args.max_cores)}, "
        f"param_grid_size={max(1, args.param_grid_size)}, "
        f"search_strategy={args.search_strategy}, "
        f"lhs_scale_mode={args.lhs_scale_mode}, "
        f"xgb_device={'auto' if 'xgb' in requested_models else 'n/a'}, "
        f"data_reader_backend={args.data_reader_backend}, "
        f"data_path={args.data_path}, "
        f"group_filter={args.only_group}, "
        f"logit_regression={args.logit_regression}, "
        f"outer_splits={args.outer_splits}, inner_splits=derived_from_outer, "
        f"calibrate_classifiers={args.calibration}, "
        f"output_level={args.output_level}, "
        f"log_level={'debug' if args.debug else 'info'}, "
        f"optuna_backend={args.optuna}, "
        f"optuna_sampler={args.optuna_sampler if args.optuna else 'n/a'}, "
        f"models={list(requested_models)}",
        level="info",
    )

    # Must be configured before importing polars for the first time.
    if args.data_reader_backend in {"auto", "polars"}:
        if "POLARS_MAX_THREADS" not in os.environ:
            filter_fraction = estimate_filter_fraction_for_path(
                args.data_path,
                filter_event_type=args.only_event_type,
                filter_transcript_filter=args.only_transcript_filter,
                filter_variability=args.only_variability,
            )
            polars_threads = recommend_polars_max_threads(
                max_cores=max(1, args.max_cores),
                filter_fraction=filter_fraction,
            )
            os.environ["POLARS_MAX_THREADS"] = str(polars_threads)
            frac_msg = (
                ""
                if filter_fraction is None
                else f", sampled_filter_fraction={filter_fraction:.3f}"
            )
            vlog(
                args.verbose,
                (
                    "Configured POLARS_MAX_THREADS="
                    f"{polars_threads} (memory-aware heuristic, max_cores={max(1, args.max_cores)}{frac_msg})"
                ),
                level="info",
            )
        else:
            vlog(
                args.verbose,
                f"Using pre-set POLARS_MAX_THREADS={os.environ['POLARS_MAX_THREADS']}",
                level="info",
            )

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
        # if smoke mode change the default to SMOKE_MODEL_TYPES, i.e., check whether requested_models are DEFAULT_MODEL_TYPES or not, if not, then use requested_models, else use SMOKE_MODEL_TYPES
        if models == DEFAULT_MODEL_TYPES:
            # extract all available models in smoke mode for comprehensive testing, regardless of user default
            models = SMOKE_MODEL_TYPES
        # Use requested models if they are not the default

    xgb_use_gpu = None
    if "xgb" in models:
        xgb_use_gpu = xgb_gpu_available()
        vlog(
            args.verbose,
            f"XGBoost runtime device check: {'cuda' if xgb_use_gpu else 'cpu'}",
            level="info",
        )

    cuml_use_gpu = None
    if "svm" in models:
        cuml_use_gpu = cuml_gpu_available()
        vlog(
            args.verbose,
            f"cuML SVM runtime device check: {'cuda' if cuml_use_gpu else 'cpu'}",
            level="info",
        )

    if "mlp" in models or "tabicl" in models:
        try:
            import torch as _torch

            _torch_device = "cpu"
            try:
                if _torch.accelerator.is_available():
                    _torch_device = str(_torch.accelerator.current_accelerator())
            except AttributeError:
                if _torch.cuda.is_available():
                    _torch_device = "cuda"
                elif (
                    hasattr(_torch.backends, "mps")
                    and _torch.backends.mps.is_available()
                ):
                    _torch_device = "mps"
        except ImportError:
            _torch_device = "torch not installed"
        if "mlp" in models:
            vlog(args.verbose, f"MLP runtime device check: {_torch_device}", level="info")
        if "tabicl" in models:
            vlog(args.verbose, f"TabICL runtime device check: {_torch_device}", level="info")

    outer_splits = args.outer_splits
    if args.smoke:
        outer_splits = 4
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
        optuna_backend=args.optuna,
        optuna_sampler=args.optuna_sampler,
        optuna_n_startup_trials=args.optuna_n_startup_trials,
        optuna_wandb_callback=args.optuna_wandb_callback,
        xgb_use_gpu=xgb_use_gpu,
        cuml_use_gpu=cuml_use_gpu,
        verbose=args.verbose,
        include_models=models,
        run_regression=not args.skip_regression,
        run_classification=not args.skip_classification,
        use_logit_regression=args.logit_regression,
        categorical_missing_strategy=args.categorical_missing_strategy,
        categorical_missing_token=args.categorical_missing_token,
        generate_html_reports=not args.no_html_reports,
        smoke_mode=args.smoke,
        smoke_max_rows=args.smoke_max_rows,
        calibrate_classifiers=args.calibration,
        tune_threshold=args.tune_threshold,
        output_level=args.output_level,
        log_level=("debug" if args.debug else "info"),
        use_wandb=args.wandb,
        wandb_project=args.wandb_project,
        wandb_entity=args.wandb_entity,
        wandb_require_auth=not args.wandb_no_require_auth,
        wandb_log_fold_table=not args.wandb_no_fold_table,
        wandb_log_tuning_details=not args.wandb_no_tuning_details,
        wandb_log_baseline_metrics=not args.wandb_no_baseline_metrics,
        wandb_fold_subruns=args.wandb_fold_subruns,
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
        f"log_level={run_cfg.log_level}, "
        f"optuna_backend={run_cfg.optuna_backend}, "
        f"optuna_sampler={run_cfg.optuna_sampler if run_cfg.optuna_backend else 'n/a'}, "
        f"models={list(run_cfg.include_models)}",
        level="info",
    )
    if run_cfg.use_wandb:
        vlog(
            args.verbose,
            "W&B settings: "
            f"project={run_cfg.wandb_project}, "
            f"entity={run_cfg.wandb_entity}, "
            f"require_auth={run_cfg.wandb_require_auth}, "
            f"fold_table={run_cfg.wandb_log_fold_table}, "
            f"tuning_details={run_cfg.wandb_log_tuning_details}, "
            f"baseline_metrics={run_cfg.wandb_log_baseline_metrics}, "
            f"fold_subruns={run_cfg.wandb_fold_subruns}",
            level="info",
        )

    # Thread limits are enforced by SLURM cgroups (-c N); no explicit guardrails needed.
    # UNUSED: Legacy explicit thread guardrail wiring kept as reference.
    # TODO: remove this dead block if cgroup-only scheduling remains the project default.
    # thread_env = _configure_parallelism_guardrails(run_cfg.max_cores)
    # vlog(args.verbose, f"Parallelism guardrails: {thread_env}", level="info")

    bundle = run_all_configurations(run_cfg)
    vlog(
        args.verbose,
        f"Completed {len(bundle['results'])} config-task runs",
        level="info",
    )


if __name__ == "__main__":
    main()
