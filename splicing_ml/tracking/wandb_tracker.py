from __future__ import annotations

"""W&B tracking adapter with a NullTracker fallback.

This module keeps experiment tracking optional and non-fatal. Call sites can
use a single tracker object without guard conditionals.
"""

import contextlib
import json
import os
import warnings
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


__all__ = ["NullTracker", "WandbTracker", "make_tracker"]


class NullTracker:
    """No-op tracker matching the W&B tracker interface."""

    def start_run(self, *args: Any, **kwargs: Any) -> None:
        """Start run."""
        return None

    def start_config_task(self, *args: Any, **kwargs: Any) -> None:
        """Start config task."""
        return None

    def log_fold_result(self, *args: Any, **kwargs: Any) -> None:
        """Log fold result."""
        return None

    def log_config_task_summary(self, *args: Any, **kwargs: Any) -> None:
        """Log config task summary."""
        return None

    def log_artifact(self, *args: Any, **kwargs: Any) -> None:
        """Log artifact."""
        return None

    def finish_config_task(self, *args: Any, **kwargs: Any) -> None:
        """Finish config task."""
        return None

    def finish_run(self, *args: Any, **kwargs: Any) -> None:
        """Finish run."""
        return None


class WandbTracker:
    """W&B-backed tracker.

    Creates one parent orchestrator run and one child run per
    (subset configuration x task) workload.
    """

    def __init__(self, wandb_module: Any, run_cfg: Any):
        """Initialize a WandbTracker instance."""
        self._wandb = wandb_module
        self._project = str(getattr(run_cfg, "wandb_project", "splicing-ml"))
        self._entity = getattr(run_cfg, "wandb_entity", None)
        self._run_cfg = run_cfg
        self._parent_run = None
        self._child_run = None
        self._child_context: dict[str, Any] = {}
        self._fold_rows: list[dict[str, Any]] = []
        self._cv_results_by_model: dict[str, list[dict[str, Any]]] = {}
        # Pooled across outer folds (one entry per model) -- used to log ONE
        # aggregated confusion matrix / ROC / PR curve per model at
        # log_config_task_summary() time, instead of one per fold. Five
        # near-duplicate single-fold curves per model was the main source of
        # "not helpful" clutter in the old per-fold logging.
        self._oof_by_model: dict[str, dict[str, list]] = {}
        # feature_name -> list of |importance| values, one per fold, per model.
        # Aggregated into one mean-importance chart per model instead of
        # wandb.sklearn.plot_feature_importances() being called once per fold
        # (that helper only supports tree/linear coef_ models anyway, and one
        # panel per fold made cross-fold stability impossible to read).
        self._importances_by_model: dict[str, dict[str, list[float]]] = {}

    @staticmethod
    def _safe_float(value: Any) -> float | None:
        """Internal helper for safe float."""
        if not isinstance(value, (int, float)):
            return None
        fval = float(value)
        if not np.isfinite(fval):
            return None
        return fval

    def _extract_tuning_summary(self, tuning_info: dict[str, Any]) -> dict[str, Any]:
        """Internal helper for extract tuning summary."""
        best_score = self._safe_float(tuning_info.get("best_score"))
        best_params = tuning_info.get("best_params", {})
        if not isinstance(best_params, dict):
            best_params = {}
        param_keys = sorted(str(k) for k in best_params.keys())
        summary: dict[str, Any] = {
            "best_score": best_score,
            "best_param_count": int(len(best_params)),
            "best_param_keys": ",".join(param_keys),
            "best_params_json": json.dumps(best_params, sort_keys=True, default=str),
        }
        if "tree_es_applied" in tuning_info:
            summary["tree_es_applied"] = bool(tuning_info.get("tree_es_applied"))
        if "tree_es_n_estimators" in tuning_info:
            n_estimators = tuning_info.get("tree_es_n_estimators")
            if isinstance(n_estimators, (int, float)):
                summary["tree_es_n_estimators"] = int(n_estimators)
        if "beta_convergence_risk" in tuning_info:
            summary["beta_convergence_risk"] = json.dumps(
                tuning_info.get("beta_convergence_risk"),
                sort_keys=True,
                default=str,
            )
        return summary

    def _log_fold_subrun(
        self, payload: dict[str, Any], fold_row: dict[str, Any]
    ) -> None:
        """Internal helper for log fold subrun."""
        if not bool(getattr(self._run_cfg, "wandb_fold_subruns", False)):
            return
        if self._child_run is None:
            return
        child_name = str(getattr(self._child_run, "name", "subset-task"))
        run_name = (
            f"{child_name}-fold{int(fold_row['fold_id']):02d}-{fold_row['model_name']}"
        )
        subrun = None
        with contextlib.suppress(Exception):
            subrun = self._wandb.init(
                project=self._project,
                entity=self._entity,
                group=str(getattr(self._child_run, "group", None)),
                job_type="subset-task-fold",
                name=run_name,
                tags=["fold", str(fold_row.get("task", ""))],
                config={
                    "fold_id": int(fold_row["fold_id"]),
                    "model_name": str(fold_row["model_name"]),
                    "task": str(fold_row["task"]),
                    "subset_context": dict(self._child_context),
                },
                reinit=True,
            )
            subrun.log(payload)
            subrun.summary.update(
                {
                    "primary_metric": fold_row.get("primary_metric"),
                    "primary_metric_value": fold_row.get("primary_metric_value"),
                    "n_test_samples": fold_row.get("n_test_samples"),
                }
            )
            subrun.finish()

    def _log_fold_table(self) -> None:
        """Internal helper for log fold table."""
        if self._child_run is None:
            return
        if not bool(getattr(self._run_cfg, "wandb_log_fold_table", True)):
            return
        if not self._fold_rows:
            return

        columns = [
            "fold_id",
            "model_name",
            "task",
            "primary_metric",
            "primary_metric_value",
            "train_primary_metric_value",
            "generalization_gap",
            "n_test_samples",
            "n_train_samples",
            "best_inner_score",
            "best_param_count",
            "best_param_keys",
            "best_params_json",
            "baseline_primary_metric_value",
            "delta_vs_baseline_primary",
        ]
        table = self._wandb.Table(columns=columns)
        for row in self._fold_rows:
            table.add_data(
                row.get("fold_id"),
                row.get("model_name"),
                row.get("task"),
                row.get("primary_metric"),
                row.get("primary_metric_value"),
                row.get("train_primary_metric_value"),
                row.get("generalization_gap"),
                row.get("n_test_samples"),
                row.get("n_train_samples"),
                row.get("best_inner_score"),
                row.get("best_param_count"),
                row.get("best_param_keys"),
                row.get("best_params_json"),
                row.get("baseline_primary_metric_value"),
                row.get("delta_vs_baseline_primary"),
            )
        with contextlib.suppress(Exception):
            self._child_run.log({"folds/by_model": table})

    def start_run(self, run_cfg: Any) -> None:
        """Start run."""
        tags = ["splicing-ml", "orchestrator"]
        if bool(getattr(run_cfg, "smoke_mode", False)):
            tags.append("smoke")
        tags.append(f"log_level:{getattr(run_cfg, 'log_level', 'info')}")
        self._parent_run = self._wandb.init(
            project=self._project,
            entity=self._entity,
            job_type="orchestrator",
            config=asdict(run_cfg),
            tags=tags,
            reinit=True,
        )

    def _delete_existing_run(self, name: str) -> None:
        """Delete any prior FINISHED/dead run with this exact display name in
        this project.

        Rerunning the same subset config produces the same deterministic
        run_name (see start_config_task) -- without this, each rerun piles up
        as a new run instead of replacing the previous one's results.

        Skips any run still in state "running" (real risk, not theoretical:
        this pipeline gets retried/resubmitted after fixes constantly --
        found 2026-07-16 during a design review, before this ever became the
        default). SLURM cancellation isn't instant, so a retry's
        start_config_task can otherwise race an old attempt that's still
        alive and still calling self._child_run.log(...) -- deleting that
        run out from under it, silently, since every W&B call site here is
        wrapped in contextlib.suppress(Exception).
        """
        with contextlib.suppress(Exception):
            api = self._wandb.Api()
            entity = self._entity or api.default_entity
            path = f"{entity}/{self._project}"
            for run in api.runs(path, filters={"display_name": name}):
                if getattr(run, "state", None) == "running":
                    continue
                run.delete()

    def start_config_task(
        self,
        cfg: Any,
        task: str,
        run_cfg: Any,
        n_samples: int,
        prep_details: dict[str, Any],
    ) -> None:
        """Start config task."""
        group = getattr(self._parent_run, "name", None)
        # Ablation runs (RunConfig.feature_groups) were previously indistinguishable
        # from full-feature runs in the W&B UI (same run-name pattern, no tag,
        # feature_groups buried inside a nested "run_config" config blob) --
        # promote it to a top-level tag/config key/run-name suffix so ablation
        # configs are filterable and groupable like the other split axes.
        feature_groups = getattr(run_cfg, "feature_groups", None)
        feature_groups_label = "+".join(feature_groups) if feature_groups else "all"
        run_name = (
            f"{cfg.event_type}-{cfg.transcript_filter}-{cfg.variability}-"
            f"{cfg.group_col}-{task}"
        )
        if feature_groups_label != "all":
            run_name += f"-feat_{feature_groups_label}"
        tags = [
            f"event:{cfg.event_type}",
            f"tx:{cfg.transcript_filter}",
            f"var:{cfg.variability}",
            f"group:{cfg.group_col}",
            f"task:{task}",
            f"features:{feature_groups_label}",
            "subset-task",
        ]
        if bool(getattr(run_cfg, "smoke_mode", False)):
            tags.append("smoke")
        if bool(getattr(run_cfg, "wandb_fold_subruns", False)):
            tags.append("fold-subruns")
        self._delete_existing_run(run_name)
        self._child_run = self._wandb.init(
            project=self._project,
            entity=self._entity,
            group=group,
            job_type="subset-task",
            name=run_name,
            tags=tags,
            config={
                "subset_config": asdict(cfg),
                "task": task,
                "run_config": asdict(run_cfg),
                "n_samples": int(n_samples),
                "feature_groups": list(feature_groups) if feature_groups else ["all"],
            },
            reinit=True,
        )
        self._child_context = {
            "event_type": str(cfg.event_type),
            "transcript_filter": str(cfg.transcript_filter),
            "variability": str(cfg.variability),
            "group_col": str(cfg.group_col),
            "task": str(task),
            "feature_groups": feature_groups_label,
        }
        self._fold_rows = []
        self._cv_results_by_model = {}
        self._oof_by_model = {}
        self._importances_by_model = {}
        # Without this, every logged metric is plotted against W&B's internal
        # auto-incrementing _step (i.e. call order), not the fold number — since
        # each fold triggers 2-3 separate .log() calls (scores, tuning_summary,
        # sklearn diagnostics) across multiple models, the default charts end up
        # as scrambled/sawtooth lines that don't compare across models at all.
        # Declaring fold_id as the step metric makes every chart use it as the
        # x-axis instead, regardless of call order.
        with contextlib.suppress(Exception):
            self._child_run.define_metric("fold_id")
            self._child_run.define_metric("*", step_metric="fold_id")
        with contextlib.suppress(Exception):
            prep_table = self._wandb.Table(columns=["key", "value"])
            for key, value in prep_details.items():
                prep_table.add_data(str(key), str(value))
            self._child_run.log({"preprocessing/details": prep_table})

    def log_fold_result(
        self,
        fold_id: int,
        model_name: str,
        scores: dict[str, float],
        tuning_info: dict[str, Any],
        y_true: np.ndarray | list[float],
        y_pred: np.ndarray | list[float],
        estimator: Any,
        task: str,
        x_train: Any,
        train_scores: dict[str, float] | None = None,
        primary_metric: str | None = None,
        baseline_scores: dict[str, float] | None = None,
        delta_vs_baseline: dict[str, float] | None = None,
    ) -> None:
        """Log fold result."""
        if self._child_run is None:
            return

        # Prefix model metrics so multiple model families can share one run.
        payload = {
            f"{model_name}/{k}": float(v)
            for k, v in scores.items()
            if isinstance(v, (int, float)) and np.isfinite(float(v))
        }
        payload[f"{model_name}/fold_id"] = int(fold_id)
        # Unprefixed fold_id (paired with the define_metric step_metric wiring
        # in start_config_task) so every chart in this run uses fold number as
        # its x-axis instead of call order.
        payload["fold_id"] = int(fold_id)
        # Log training scores (prefixed train/) to expose overfitting gap.
        if train_scores:
            for k, v in train_scores.items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/train_{k}"] = float(v)
        if bool(getattr(self._run_cfg, "wandb_log_baseline_metrics", True)):
            for k, v in (baseline_scores or {}).items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/baseline_{k}"] = float(v)
            for k, v in (delta_vs_baseline or {}).items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/delta_vs_baseline_{k}"] = float(v)
        self._child_run.log(payload)

        if bool(getattr(self._run_cfg, "wandb_log_tuning_details", True)):
            # Previously also logged the raw tuning_info dict wholesale
            # (f"{model_name}/tuning") — nested dicts (best_params objects,
            # full cv_results arrays) don't render as a chart, just clutter as
            # unreadable panels. The flat summary below is the useful part;
            # _log_cv_result_tables() already exposes the full cv_results as a
            # proper Table.
            with contextlib.suppress(Exception):
                tuning_summary = self._extract_tuning_summary(tuning_info)
                self._child_run.log({f"{model_name}/tuning_summary": tuning_summary})

        primary_key = str(primary_metric or "")
        primary_value = (
            self._safe_float(scores.get(primary_key)) if primary_key else None
        )
        train_primary = (
            self._safe_float((train_scores or {}).get(primary_key))
            if primary_key
            else None
        )
        gap = (
            float(train_primary - primary_value)
            if train_primary is not None and primary_value is not None
            else None
        )
        tuning_summary_for_row = self._extract_tuning_summary(tuning_info)
        baseline_primary = (
            self._safe_float((baseline_scores or {}).get(primary_key))
            if primary_key
            else None
        )
        delta_primary = (
            self._safe_float((delta_vs_baseline or {}).get(primary_key))
            if primary_key
            else None
        )
        fold_row = {
            "fold_id": int(fold_id),
            "model_name": str(model_name),
            "task": str(task),
            "primary_metric": primary_key,
            "primary_metric_value": primary_value,
            "train_primary_metric_value": train_primary,
            "generalization_gap": gap,
            "n_test_samples": int(len(y_true)) if y_true is not None else None,
            "n_train_samples": int(getattr(x_train, "shape", [len(x_train)])[0]),
            "best_inner_score": tuning_summary_for_row.get("best_score"),
            "best_param_count": tuning_summary_for_row.get("best_param_count"),
            "best_param_keys": tuning_summary_for_row.get("best_param_keys"),
            "best_params_json": tuning_summary_for_row.get("best_params_json"),
            "baseline_primary_metric_value": baseline_primary,
            "delta_vs_baseline_primary": delta_primary,
        }
        self._fold_rows.append(fold_row)
        cv_results = (
            tuning_info.get("cv_results") if isinstance(tuning_info, dict) else None
        )
        if cv_results and isinstance(cv_results, dict):
            entry = {"outer_fold": int(fold_id), **cv_results}
            self._cv_results_by_model.setdefault(model_name, []).append(entry)
        self._log_fold_subrun(payload, fold_row)

        # Accumulate out-of-fold predictions for a SINGLE pooled confusion
        # matrix / ROC / PR curve per model, logged once in
        # log_config_task_summary() -- five near-duplicate single-fold curves
        # per model was clutter, not signal (each fold is a genuinely
        # different chromosome/ontology group, so pooling them into one
        # curve is also the statistically correct thing to do here, not just
        # a display simplification).
        if task == "classification":
            oof = self._oof_by_model.setdefault(
                model_name, {"y_true": [], "y_pred": []}
            )
            oof["y_true"].extend(np.asarray(y_true, dtype=int).tolist())
            oof["y_pred"].extend(np.asarray(y_pred, dtype=float).tolist())

        # Accumulate |importance|/|coef| per feature for a SINGLE mean-across-folds
        # chart per model (logged in log_config_task_summary()), instead of
        # wandb.sklearn.plot_feature_importances() firing once per fold (only
        # supports tree/linear models anyway, and per-fold panels can't show
        # cross-fold stability).
        with contextlib.suppress(Exception):
            model = estimator.named_steps.get("model", estimator)
            prep = estimator.named_steps.get("prep")
            feature_names = prep.get_feature_names_out() if prep is not None else None
            importances = None
            if hasattr(model, "feature_importances_"):
                importances = np.asarray(model.feature_importances_, dtype=float)
            elif hasattr(model, "coef_"):
                coef = np.asarray(model.coef_, dtype=float)
                importances = np.abs(coef[0] if coef.ndim > 1 else coef)
            if importances is not None and feature_names is not None:
                bucket = self._importances_by_model.setdefault(model_name, {})
                for name, val in zip(feature_names, importances):
                    bucket.setdefault(str(name), []).append(float(val))

    def _log_model_comparison(self, primary_metric: str) -> None:
        """Log outer-fold model comparison: bar chart (mean) + fold-detail table."""
        if self._child_run is None or not self._fold_rows:
            return
        from collections import defaultdict

        scores_by_model: dict[str, list[float]] = defaultdict(list)
        for row in self._fold_rows:
            v = row.get("primary_metric_value")
            if v is not None:
                scores_by_model[str(row["model_name"])].append(float(v))
        if not scores_by_model:
            return

        summary_table = self._wandb.Table(columns=["model", "mean", "std", "n_folds"])
        for model_name, vals in sorted(scores_by_model.items()):
            summary_table.add_data(
                model_name,
                float(np.mean(vals)),
                float(np.std(vals)),
                len(vals),
            )

        detail_table = self._wandb.Table(columns=["model", "fold_id", "score"])
        for row in self._fold_rows:
            v = row.get("primary_metric_value")
            if v is not None:
                detail_table.add_data(
                    str(row["model_name"]),
                    int(row["fold_id"]),
                    float(v),
                )

        # One multi-line chart (all models, x=fold_id) instead of relying on
        # the per-model auto-generated single-line panels from log_fold_result
        # — those land at different points on W&B's default x-axis whenever
        # models are logged in different .log() calls, so they were never
        # visually comparable without this. Sort each model's (fold_id, score)
        # pairs by fold_id first so the line doesn't zigzag if folds completed
        # out of order (parallel fitting).
        fold_score_pairs_by_model: dict[str, list[tuple[int, float]]] = defaultdict(
            list
        )
        for row in self._fold_rows:
            v = row.get("primary_metric_value")
            if v is not None:
                fold_score_pairs_by_model[str(row["model_name"])].append(
                    (int(row["fold_id"]), float(v))
                )
        model_names_sorted = sorted(fold_score_pairs_by_model)
        line_chart = None
        with contextlib.suppress(Exception):
            xs, ys = [], []
            for m in model_names_sorted:
                pairs = sorted(fold_score_pairs_by_model[m])
                xs.append([p[0] for p in pairs])
                ys.append([p[1] for p in pairs])
            line_chart = self._wandb.plot.line_series(
                xs=xs,
                ys=ys,
                keys=model_names_sorted,
                title=f"{primary_metric} by fold, per model",
                xname="fold_id",
            )

        payload = {
            "comparison/model_summary": summary_table,
            "comparison/model_fold_detail": detail_table,
            "comparison/bar_chart": self._wandb.plot.bar(
                summary_table,
                "model",
                "mean",
                title=f"Model comparison — mean {primary_metric} (outer folds)",
            ),
        }
        if line_chart is not None:
            payload["comparison/line_chart"] = line_chart
        with contextlib.suppress(Exception):
            self._child_run.log(payload)

    def _log_pooled_classification_diagnostics(self) -> None:
        """Log ONE confusion matrix / ROC / PR curve per model, pooled across
        all outer folds -- replaces the old per-fold versions (5 near-duplicate
        panels per model was clutter; pooling is also the statistically correct
        choice since each fold is a genuinely different chromosome/ontology
        group, not a repeat of the same evaluation)."""
        if self._child_run is None or not self._oof_by_model:
            return
        for model_name, oof in self._oof_by_model.items():
            if not oof["y_true"]:
                continue
            with contextlib.suppress(Exception):
                yt = np.asarray(oof["y_true"], dtype=int)
                yp = np.asarray(oof["y_pred"], dtype=float)
                if len(np.unique(yt)) < 2:
                    continue
                y_probas = np.column_stack([1.0 - yp, yp])
                labels = ["0", "1"]
                hard = (yp >= 0.5).astype(int)
                self._child_run.log(
                    {
                        f"{model_name}/confusion_matrix_pooled": self._wandb.plot.confusion_matrix(
                            y_true=yt.tolist(),
                            preds=hard.tolist(),
                            class_names=labels,
                        ),
                        f"{model_name}/roc_pooled": self._wandb.plot.roc_curve(
                            yt, y_probas, labels
                        ),
                        f"{model_name}/pr_curve_pooled": self._wandb.plot.pr_curve(
                            yt, y_probas, labels
                        ),
                    }
                )

    def _log_aggregated_feature_importance(self, top_n: int = 30) -> None:
        """Log one mean-across-folds |importance|/|coef| chart per model,
        instead of wandb.sklearn.plot_feature_importances() firing once per
        fold (tree/linear models only, no cross-fold stability visible)."""
        if self._child_run is None or not self._importances_by_model:
            return
        for model_name, feature_vals in self._importances_by_model.items():
            with contextlib.suppress(Exception):
                rows = [
                    (name, float(np.mean(vals)), float(np.std(vals)), len(vals))
                    for name, vals in feature_vals.items()
                ]
                rows.sort(key=lambda r: r[1], reverse=True)
                rows = rows[:top_n]
                table = self._wandb.Table(
                    columns=["feature", "mean_importance", "std_importance", "n_folds"],
                    data=[list(r) for r in rows],
                )
                self._child_run.log(
                    {
                        f"{model_name}/feature_importance": table,
                        f"{model_name}/feature_importance_bar": self._wandb.plot.bar(
                            table,
                            "feature",
                            "mean_importance",
                            title=f"{model_name} — top {top_n} features (mean |importance|/|coef| across folds)",
                        ),
                    }
                )

    def _log_cv_result_tables(self) -> None:
        """Log one aggregated inner-CV results Table per model (across all outer folds)."""
        if self._child_run is None or not self._cv_results_by_model:
            return
        for model_name, fold_entries in self._cv_results_by_model.items():
            with contextlib.suppress(Exception):
                # Determine split column names from first entry.
                split_keys: list[str] = []
                for entry in fold_entries:
                    split_keys = entry.get("split_keys", [])
                    if split_keys:
                        break
                columns = [
                    "outer_fold",
                    "rank",
                    "mean_score",
                    "std_score",
                    "params",
                ] + split_keys
                table = self._wandb.Table(columns=columns)
                for entry in fold_entries:
                    outer_fold = entry["outer_fold"]
                    ranks = entry.get("rank_test_score", [])
                    means = entry.get("mean_test_score", [])
                    stds = entry.get("std_test_score", [])
                    params = entry.get("params_json", [])
                    split_scores = entry.get("split_scores", {})
                    n = entry.get("n_candidates", len(means))
                    for i in range(n):
                        split_vals = [
                            split_scores.get(k, [None] * n)[i] for k in split_keys
                        ]
                        table.add_data(
                            outer_fold,
                            ranks[i] if i < len(ranks) else None,
                            means[i] if i < len(means) else None,
                            stds[i] if i < len(stds) else None,
                            params[i] if i < len(params) else None,
                            *split_vals,
                        )
                self._child_run.log({f"{model_name}/cv_results": table})

    def log_config_task_summary(self, result: dict[str, Any]) -> None:
        """Log config task summary."""
        if self._child_run is None:
            return
        with contextlib.suppress(Exception):
            self._child_run.summary.update(
                {
                    "status": result.get("status"),
                    "task": result.get("task"),
                    "primary_metric": result.get("primary_metric"),
                    "primary_metric_mean": result.get("primary_metric_mean"),
                    "primary_metric_std": result.get("primary_metric_std"),
                    "primary_metric_bootstrap_ci": result.get(
                        "primary_metric_bootstrap_ci"
                    ),
                    # Per-model breakdown -- primary_metric_mean/std above blend
                    # ALL models together, which is not a meaningful number for
                    # comparing model architectures. See pipeline.py's
                    # primary_metric_by_model.
                    "primary_metric_by_model": result.get("primary_metric_by_model"),
                    "n_samples": result.get("n_samples"),
                    "n_fold_results": len(result.get("fold_results", [])),
                }
            )
        with contextlib.suppress(Exception):
            self._log_fold_table()
        with contextlib.suppress(Exception):
            self._log_model_comparison(str(result.get("primary_metric", "")))
        with contextlib.suppress(Exception):
            self._log_cv_result_tables()
        with contextlib.suppress(Exception):
            self._log_pooled_classification_diagnostics()
        with contextlib.suppress(Exception):
            self._log_aggregated_feature_importance()

    def log_artifact(self, path: str | Path, name: str, artifact_type: str) -> None:
        """Log artifact."""
        if self._parent_run is None:
            return
        p = Path(path)
        if not p.exists():
            return
        with contextlib.suppress(Exception):
            artifact = self._wandb.Artifact(name=name, type=artifact_type)
            artifact.add_file(str(p))
            self._parent_run.log_artifact(artifact)

    def finish_config_task(self) -> None:
        """Finish config task."""
        if self._child_run is None:
            return
        with contextlib.suppress(Exception):
            self._child_run.finish()
        self._child_run = None
        self._fold_rows = []
        self._cv_results_by_model = {}
        self._oof_by_model = {}
        self._importances_by_model = {}
        self._child_context = {}

    def finish_run(self, all_results: list[dict[str, Any]]) -> None:
        """Finish run.

        Master cross-run comparison table: ONE ROW PER (config, task, model),
        not per config-task. The previous version had one row per config-task
        with a single blended `primary_metric_mean` averaged across every
        model AND fold together -- meaningless for "which model wins", and
        with no `feature_groups`/`model_name` columns at all there was no way
        to filter/group by split axis, ablation, or model in this table.
        Every axis relevant to "compare models/splits/ablations" is now its
        own column, so W&B's native Table groupby/bar/box view (interactive,
        in the UI) can build any of those comparisons directly from this one
        table instead of needing a bespoke pre-baked chart per combination.
        """
        if self._parent_run is None:
            return
        with contextlib.suppress(Exception):
            summary_table = self._wandb.Table(
                columns=[
                    "event_type",
                    "transcript_filter",
                    "variability",
                    "group_col",
                    "feature_groups",
                    "task",
                    "status",
                    "model_name",
                    "primary_metric",
                    "mean",
                    "std",
                    "n_folds",
                ]
            )
            for row in all_results:
                cfg = row.get("config", {})
                base = (
                    str(cfg.get("event_type", "")),
                    str(cfg.get("transcript_filter", "")),
                    str(cfg.get("variability", "")),
                    str(cfg.get("group_col", "")),
                    "+".join(row.get("feature_groups") or ["all"]),
                    str(row.get("task", "")),
                    str(row.get("status", "")),
                )
                by_model = row.get("primary_metric_by_model") or {}
                if by_model:
                    for model_name, stats in by_model.items():
                        summary_table.add_data(
                            *base,
                            str(model_name),
                            str(row.get("primary_metric", "")),
                            float(stats.get("mean", np.nan)),
                            float(stats.get("std", np.nan)),
                            int(stats.get("n_folds", 0)),
                        )
                else:
                    # failed/skipped config-tasks have no per-model breakdown --
                    # keep one row so they're still visible in the table.
                    summary_table.add_data(
                        *base,
                        "",
                        str(row.get("primary_metric", "")),
                        float(row.get("primary_metric_mean", np.nan)),
                        float(row.get("primary_metric_std", np.nan)),
                        0,
                    )
            self._parent_run.log({"summary/by_config_task": summary_table})

        # Three quick rankings answering "which model/split/ablation wins",
        # faceted by event_type (found 2026-07-16 during a design review,
        # before this ever became the default: RI vs SE have such different
        # baseline difficulty/scale -- demonstrated repeatedly this session
        # -- that a single mean blended across both event types is the exact
        # same anti-pattern primary_metric_by_model was introduced to fix for
        # the per-model blend, just recurring one level up. `summary/by_config_task`
        # above is the rich table for anything even more specific (filter to
        # one transcript_filter, one variability, etc.) via W&B's own
        # interactive Table UI.
        event_types = sorted(
            {str(row.get("config", {}).get("event_type", "")) for row in all_results}
            - {""}
        )
        for group_col_name, chart_title in (
            ("model_name", "Mean primary metric by MODEL"),
            ("feature_groups", "Mean primary metric by FEATURE GROUP (ablations)"),
            ("group_col", "Mean primary metric by CV SPLIT AXIS (seqnames vs ontology)"),
        ):
            for event_type in event_types:
                with contextlib.suppress(Exception):
                    buckets: dict[str, list[float]] = {}
                    for row in all_results:
                        cfg = row.get("config", {})
                        if str(cfg.get("event_type", "")) != event_type:
                            continue
                        by_model = row.get("primary_metric_by_model") or {}
                        for model_name, stats in by_model.items():
                            key = (
                                model_name
                                if group_col_name == "model_name"
                                else (
                                    "+".join(row.get("feature_groups") or ["all"])
                                    if group_col_name == "feature_groups"
                                    else str(cfg.get("group_col", ""))
                                )
                            )
                            mean_val = stats.get("mean")
                            if mean_val is not None and np.isfinite(mean_val):
                                buckets.setdefault(str(key), []).append(float(mean_val))
                    if not buckets:
                        continue
                    rank_table = self._wandb.Table(
                        columns=[group_col_name, "mean", "n"]
                    )
                    for key, vals in sorted(buckets.items()):
                        rank_table.add_data(key, float(np.mean(vals)), len(vals))
                    self._parent_run.log(
                        {
                            f"ranking/by_{group_col_name}__{event_type}": rank_table,
                            f"ranking/by_{group_col_name}__{event_type}_bar": self._wandb.plot.bar(
                                rank_table,
                                group_col_name,
                                "mean",
                                title=f"{chart_title} ({event_type})",
                            ),
                        }
                    )

        with contextlib.suppress(Exception):
            self._parent_run.finish()
        self._parent_run = None


def make_tracker(run_cfg: Any) -> NullTracker | WandbTracker:
    """Create a tracker from runtime config with graceful fallback."""
    if not bool(getattr(run_cfg, "use_wandb", False)):
        return NullTracker()

    try:
        import wandb
    except Exception:
        warnings.warn(
            "use_wandb=True but wandb is unavailable; proceeding with NullTracker",
            RuntimeWarning,
            stacklevel=2,
        )
        return NullTracker()

    project = str(getattr(run_cfg, "wandb_project", "splicing-ml")).strip()
    if not project:
        raise ValueError("W&B project must be non-empty when --wandb is enabled")

    require_auth = bool(getattr(run_cfg, "wandb_require_auth", True))
    has_env_api_key = bool(os.environ.get("WANDB_API_KEY"))
    has_stored_auth = False
    if not has_env_api_key:
        # Accept credentials loaded via `wandb login` (e.g., ~/.netrc).
        with contextlib.suppress(Exception):
            has_stored_auth = bool(getattr(wandb.Api(), "api_key", None))

    has_auth = has_env_api_key or has_stored_auth
    if require_auth and not has_auth:
        raise RuntimeError(
            "W&B is enabled but no credentials were found. "
            "Run `wandb login`, export WANDB_API_KEY, or disable strict auth "
            "with --wandb-no-require-auth."
        )
    if not require_auth and not has_auth:
        warnings.warn(
            "use_wandb=True but W&B credentials are missing; proceeding with NullTracker "
            "because strict auth is disabled",
            RuntimeWarning,
            stacklevel=2,
        )
        return NullTracker()

    return WandbTracker(wandb_module=wandb, run_cfg=run_cfg)
